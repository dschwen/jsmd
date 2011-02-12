var jsmd = (function(){
  // Vector constructor
  function Vector(x,y) {
    this.x = x || 0.0;
    this.y = y || 0.0;
  }
  Vector.prototype.len = function() {
    return Math.sqrt( this.x * this.x + this.y * this.y );
  }
  Vector.prototype.len2 = function() {
    return this.x * this.x + this.y * this.y;
  }
  Vector.prototype.pbclen = function(w,h) {
    var dx = this.x - Math.round(this.x/w) * w,
        dy = this.y - Math.round(this.y/h) * h;
    return Math.sqrt( dx*dx + dy*dy );
  }
  Vector.prototype.normalize = function() {
    var len = this.len();
    if( len === 0 ) {
      this.x = 1; this.y = 0;
    } else {
      this.x /= len; this.y /= len;
    }
  }
  Vector.distance = function(a,b) {
    return Math.sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) );
  }
  Vector.distance2 = function(a,b) {
    return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
  }
  Vector.pbcdistance = function(a,b,w,h) {
    var dx = a.x-b.x,
        dy = a.y-b.y;
    dx -= Math.round(dx/w) * w;
    dy -= Math.round(dy/h) * h;
    return Math.sqrt( dx*dx + dy*dy );
  }
  Vector.pbcdistance2 = function(a,b,w,h) {
    var dx = a.x-b.x,
        dy = a.y-b.y;
    dx -= Math.round(dx/w) * w;
    dy -= Math.round(dy/h) * h;
    return dx*dx + dy*dy;
  }
  Vector.prototype.wrap = function(w,h) {
    this.x -= Math.floor(this.x/w) * w;
    this.y -= Math.floor(this.y/h) * h;
  }
  Vector.prototype.dwrap = function(w,h) {
    this.x -= Math.round(this.x/w) * w;
    this.y -= Math.round(this.y/h) * h;
  }
  Vector.sub = function(a,b) {
    return new Vector( a.x-b.x, a.y-b.y );
  }
  Vector.prototype.sub = function(a) {
    this.x -= a.x;
    this.y -= a.y;
  }
  Vector.add = function(a,b) {
    return new Vector( a.x+b.x, a.y+b.y );
  }
  Vector.prototype.add = function(a) {
    this.x += a.x;
    this.y += a.y;
  }
  Vector.scale = function(a,b) {
    return new Vector( a.x*b, a.y*b );
  }
  Vector.prototype.scale = function(a) {
    this.x *= a;
    this.y *= a;
  }
  Vector.smul = function(a,b) { // scalar product
    return a.x*b.x + a.y*b.y;
  }
  Vector.prototype.smul = function(a) { // scalar product
    return this.x*a.x + this.y*a.y;
  }
  Vector.prototype.proj = function(a) { // assume a is a unit vector!
    var s = this.smul(a);
    this.x = a.x * s;
    this.y = a.y * s;
  }
  Vector.prototype.zero = function() {
    this.x = 0.0;
    this.y = 0.0;
  }
  Vector.prototype.set = function(a) {
    this.x = a.x;
    this.y = a.y;
  }
  Vector.random = function(w,h) {
    return new Vector( Math.random()*w, Math.random()*h );
  }

  // Atom constructor
  function Atom(x,y) {
    this.p = new Vector(x,y);
    this.v = new Vector(0,0);
    this.f = new Vector(0,0);
    this.t = 0;
  }

  // AtomType constructor
  function AtomType() {
    this.r = 3;
    this.color  = "rgb(255,0,0)";
    this.Z = 1;
    this.m = 1.0;
  }

  // Barrier constructor
  function Barrier(x1,y1,x2,y2) {
    this.p = [ new Vector(x1,y1),  new Vector(x2,y2) ];

    // normal vector
    this.n = new Vector( y2-y1, x1-x2 );
    this.n.normalize();

    // colinear vector
    this.c = new Vector( x2-x1, y2-y1 );
    this.l = this.c.pbclen(); // length
    this.c.scale(1.0/this.l);

    // default interation type
    this.t = 0;
  }
  Barrier.prototype.dist = function(a) {
    var sv = new Vector( this.p[0].x - a.x, this.p[0].y - a.y ),
        s = -(sv.x*this.c.x + sv.y*this.c.y);

    if( s <= 0 ) { // return vector to endpoint 0
      return sv;
    } else if( s >= this.l ) { // return vector to endpoint 1
      return new Vector( this.p[1].x - a.x, this.p[1].y - a.y );
    } else { // return vector to line
      sv.proj( this.n );
      return sv;
    }
  }

  // Simulation constructor (this contains all the important logic)
  function Simulation(w,h, options) {
    this.atoms    = []; // list of atoms
    this.barriers = []; // list of barriers
    this.types    = []; // list of atom types

    // interaction matrix
    this.interaction = [];

    // box dimensions
    this.w = w;
    this.h = h;

    // canvas for visualization
    this.canvas = {};

    // setup default render chain
    this.renderChain = [ renderAtoms, renderForces, renderBarriers ];

    // linkcell data
    this.lc = {};

    // neighborlist data
    this.nl = { dr : 0.0, data : [], count : 0 };

    // timestep data
    this.dt = 0.0;   // set by dynamic timestepper
    this.step = 0;   // number of current MD step
    this.time = 0.0; // expired simulation time
  }
  Simulation.prototype.setInteraction = function(t,f) {
    // set the intercation function f(r,t) for the t=[a,b] atom types
    for( var i = 0; i < 2; ++i ) {
      if( this.interaction[t[i]] === undefined ) {
        this.interaction[t[i]] = [];
      }
      this.interaction[t[i]][t[1-i]] = f;
    }
  }
  Simulation.prototype.setCutOff = function(rc,rp) {
    // initialize linkcells
    var rm = rc+rp,
        i, l;
    this.rm2 = rm*rm; // safe distance squared (for fast neighborlist builds)
    this.rp = rp;
    this.rc = rc;
    this.rc2 = rc*rc;

    // number and size of cells (have at least 3 cells in each direction)
    this.lc.nx = Math.max( 3, Math.floor(this.w/rm) );
    this.lc.ny = Math.max( 3, Math.floor(this.h/rm) );
    this.lc.dx = this.w/this.lc.nx;
    this.lc.dy = this.h/this.lc.ny;

    // build the 2D linkcell grid
    l = new Array(this.lc.nx);
    for( i = 0; i < this.lc.nx; ++i ) {
      l[i] = new Array(this.lc.ny);
    }
    this.lc.data = l;
  }
  Simulation.prototype.clearLinkcell = function() {
    // clear each cell
    var i,j;
    for( i = 0; i < this.lc.nx; ++i ) {
      for( j = 0; j < this.lc.ny; ++j ) {
        this.lc.data[i][j] = [];
      }
    }
  }
  Simulation.prototype.updateLinkcell = function() {
    // clear first
    this.clearLinkcell();

    // repopulate with all atoms
    var lx, ly, i;
    for( i = 0; i < this.atoms.length; ++i ) {
      lx = Math.floor(this.atoms[i].p.x/this.lc.dx) % this.lc.nx;
      ly = Math.floor(this.atoms[i].p.y/this.lc.dy) % this.lc.ny;
      this.lc.data[lx][ly].push(i);
    }
  }
  Simulation.prototype.clearNeighborlist = function() {
    for( var i = 0; i < this.atoms.length; ++i ) {
      this.nl.data[i] = [];
    }
  }
  Simulation.prototype.updateNeighborlist = function(dr) {
    // check if update is necessary (call without parameter to force update)
    if( dr !== undefined ) {
      this.nl.dr += 2.0*dr;
      if( this.nl.dr > this.rp ) { 
        return;
      }
    }

    // update Linkcells first
    this.updateLinkcell()

    // clear Neigborlist
    this.clearNeighborlist()

    var i,j,i2,j2,k,l,m,
        ka, la, ll,
        n = [ [0,1],[1,this.lc.ny-1],[1,0],[1,1] ]; // half the neighbors

    // loop over all cells
    for( i = 0; i < this.lc.nx; ++i ) {
      for( j = 0; j < this.lc.ny; ++j ) {
        // loop over all atoms in local cell
        ll = this.lc.data[i][j].length;
        for( k = 0; k < ll; ++k ) {
          // loop over remaining atoms in local cell
          ka = this.lc.data[i][j][k];
          for( l = k+1; l < ll; ++l ) {
            la = this.lc.data[i][j][l];
            if( jsmd.Vector.pbcdistance2( this.atoms[ka].p, this.atoms[la].p, this.w, this.h ) < this.rm2 ) {
              this.nl.data[ka].push(la);
            }
          }

          // loop over half the neighbor cells
          for( m = 0; m < n.length; ++m ) {
            i2 = (i+n[m][0]) % this.lc.nx;
            j2 = (j+n[m][1]) % this.lc.ny;
            // loop over all atoms in that neighbor cells
            for( l = 0; l < this.lc.data[i2][j2].length; ++l ) {
              la = this.lc.data[i2][j2][l];
              if( jsmd.Vector.pbcdistance2( this.atoms[ka].p, this.atoms[la].p, this.w, this.h ) < this.rm2 ) {
                this.nl.data[ka].push(la);
              }
            }
          }
        }
        // end local cell
      }
    }

    // increase update counter
    this.nl.count++;
  }
  Simulation.prototype.updateForces = function() {
    var i,j,k,  // integer
        F,dr, // float
        f, // function
        rvec; // Vector

    // zero forces set timestep
    for( i = 0; i < this.atoms.length; ++i ) {
      this.atoms[i].f.x = 0.0;
      this.atoms[i].f.y = 0.0;
    }

    // build up forces (newtonian)
    for( i = 0; i < this.atoms.length; ++i ) {
      //for( j = i+1; j < atoms.length; ++j ) {     // old iteration over all other atoms
      for( k = 0; k < this.nl.data[i].length; ++k ) {  // new iteration over neighborlist
        j = this.nl.data[i][k];

        rvec = jsmd.Vector.sub( this.atoms[j].p, this.atoms[i].p);
        rvec.dwrap( this.w, this.h );
        dr = rvec.len2();
        if( dr < this.rc2 ) {
          dr = Math.sqrt(dr);
          f = this.interaction[this.atoms[i].t][this.atoms[j].t];
          if( f !== undefined ) {
            F = f.call( this, dr, [ this.atoms[i].t, this.atoms[j].t ] );
            rvec.scale(F/dr);
            this.atoms[i].f.add(rvec);
            this.atoms[j].f.sub(rvec);
          }
        }
      }
    }

    // add barrier interaction
    for( i = 0; i < this.atoms.length; ++i ) {
      for( j = 0; j < this.barriers.length; ++j ) {
        // find distance to barrier
        rvec = this.barriers[j].dist(this.atoms[i].p);
        dr = rvec.pbclen();
        if( dr < this.rc ) {
          f = this.interaction[this.atoms[i].t][this.barriers[j].t];
          if( f !== undefined ) {
            F = f.call( this, dr, [ this.atoms[i].t, this.barriers[j].t ] );
            //$('#mdlog').text(dr+','+f);
            rvec.scale(F/dr);
            this.atoms[i].f.add(rvec);
          }
        }
      }
    }
  }
  Simulation.prototype.velocityVerlet = function() {
    var i,
        dt = this.dt, m,
        dmax, rmax = 0.0, vmax = 0.0, amax = 0.0,
        dp = new jsmd.Vector();

    // first velocity verlet step
    for( i = 0; i < this.atoms.length; ++i ) {
      m = this.types[this.atoms[i].t].m;
      dp.set( this.atoms[i].v ); dp.scale(dt); // dp = v*dt
      dp.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*dt*dt) );
      this.atoms[i].p.add(dp);
      this.atoms[i].p.wrap( this.w, this.h );
      this.atoms[i].v.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*dt) );

      // track maximum displacement
      rmax = Math.max( rmax, dp.len2() );
    }

    this.updateNeighborlist(Math.sqrt(rmax));
    this.updateForces();

    // second velocity verlet step
    for( i = 0; i < this.atoms.length; ++i ) {
      this.atoms[i].v.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*dt) );

      // linear drag term
      this.atoms[i].v.scale(0.995);

      // maximum velocity and acceleration
      vmax = Math.max( vmax, this.atoms[i].v.len2() );
      amax = Math.max( amax, this.atoms[i].f.len2()/(0.25*m*m) );
    }

    // increase step counters
    this.step++;
    this.time += dt;

    // compute timestep
    dmax = 0.025;
    vmax = Math.sqrt(vmax);
    amax = Math.sqrt(amax);
    this.dt = Math.min( 0.01, dmax/vmax, Math.sqrt(2*dmax/amax) );
  }
  Simulation.prototype.setCanvas = function(canvas) {
    this.canvas = {
      node : canvas,
      w : canvas.width,
      h : canvas.height,
      ctx : canvas.getContext('2d')
    };
  }
  Simulation.prototype.draw = function() {
    var c = this.canvas.ctx,
        i;

    // setup transform
    c.setTransform( this.canvas.w/this.w, 0,0, this.canvas.h/this.h, 0,0 );

    // clear
    c.fillStyle = "rgb(200,200,200)";
    c.fillRect(0, 0, 800, 500);

    // process render chain
    for( i = 0; i < this.renderChain.length; ++i ) {
      this.renderChain[i].call(this,c);
    }
  }


  //
  // built-in render routines
  //

  // render atoms as disks with radii and colors determined by type
  function renderAtoms(c) {
    var x,y,r,i;

    function drawAtom(x,y,r) {
      c.beginPath();
      c.arc( x, y, r, 0, Math.PI*2.0, true);
      c.closePath();
      c.fill();
    }

    // draw nuclei
    c.strokeStyle = "rgba(0,100,0,0.5)";
    for( i = 0; i < this.atoms.length; ++i ) {
      c.fillStyle = this.types[this.atoms[i].t].color;
      x = this.atoms[i].p.x;
      y = this.atoms[i].p.y;
      r = this.types[this.atoms[i].t].r;

      // draw atom
      drawAtom(x,y,r);

      // draw wrap-around copies
      if( x <= r || y <= r || x+r > this.w || y+r > this.h ) {
        if( x <= r ) drawAtom(x+this.w,y,r);
        if( y <= r ) drawAtom(x,y+this.h,r);
        if( x <= r &&  y <= r ) drawAtom(x+this.w,y+this.h,r);
        if( x+r > this.w ) drawAtom(x-this.w,y,r);
        if( y+r > this.h ) drawAtom(x,y-this.h,r);
        if( x+r > this.w && y+r > this.h) drawAtom(x-this.w,y-this.h,r);
      }
    }
  }

  // draw barriers as lines
  function renderBarriers(c) {
    // draw barriers
    var i;
    c.strokeStyle = "rgba(0,0,100,1)";
    for( i = 0; i < this.barriers.length; ++i ) {
      c.beginPath();
      c.moveTo( this.barriers[i].p[0].x, this.barriers[i].p[0].y );
      c.lineTo( this.barriers[i].p[1].x, this.barriers[i].p[1].y );
      c.closePath();
      c.stroke();
    }
  }

  // visualize forces acting on individual atoms
  function renderForces(c) {
    // draw forces
    var i;
    c.strokeStyle = "rgba(0,100,0,0.5)";
    c.lineWidth = 0.01;
    for( i = 0; i < this.atoms.length; ++i ) {
      c.beginPath();
      c.moveTo( this.atoms[i].p.x, this.atoms[i].p.y );
      c.lineTo( this.atoms[i].p.x + this.atoms[i].f.x * 0.01, this.atoms[i].p.y + this.atoms[i].f.y * 0.01 );
      c.closePath();
      c.stroke();
    }
  }


  //
  // built-in force functions
  //

  // Lennard-Jones potential for pair equilibrium distance re and well depth e
  function forceLJ( rm_, e_ ) {
    var e  = e_  !== undefined ? e_ : 10.0,
        rm = rm_ !== undefined ? rm_ : 1.0,
        rm6 = Math.pow(rm,6),
        rm12 = rm6*rm6;

    function calculateForce(r) {
      r /= 4;
      return 12.0*e*( rm6*Math.pow(r,-7.0) - rm12*Math.pow(r,-13.0) );
    }

    return calculateForce;
  }

  // Morse potential for pair equilibrium distance re, potential softness 1/a, and well depth De
  function forceMorse( re_, a_, De_ ) {
    //De*( 1-exp(-a*(x-re)) )**2, (x-re)**2
    var re = re_ !== undefined ? re_ : 1.5,
        a  = a_  !== undefined ? a_  : 2.0,
        De = De_ !== undefined ? De_ : 1.0;

    function calculateForce(r) {
      var ex = Math.exp( -a*(r-re) );
      return 2.0 * a * De * (1-ex) * ex;
    }

    return calculateForce;
  }

  // return a force linearly interpolated from tabulated values of function f (fast!)
  function forceTabulated(f, dr_, rc_ ) {
    // parameters and pretabulated values stored in closure
    var dr = dr_ !== undefined ? dr_ : 0.01,
        rc = rc_ !== undefined ? rc_ : 10.0,
        mul = rc/dr,
        table = [],
        i;
    for( i = 0; i <= mul*rc+10; ++i ) {
      table[i] = f.call( this, i/mul );
    }

    // interpolation function (called from updateForces)
    function interpolate(r,t) {
      r *= mul;
      var b = Math.floor(r);
      r -= b;
      return (1-r)*table[b] + r*table[b+1];
    }

    return interpolate;
  }

  // return numerical derivative of f
  function forceNumericalDiff( f, dr_ ) {
    var dr = dr_ !== undefined ? dr_ : 0.0001;
    
    function differentiate(r) {
       return 0.5 * ( f.call(this,r-dr) + f.call(this,r+dr) ) / dr;
    }

    return differentiate;
  }

  // ZBL fore function for nuclear charges Z1 and Z2
  function forceZBL( Z1_, Z2_ ) {
  }

  // return force splined together from f1 for r<=r1, a spline for r1<r<=r2, and f2 for r2<r
  function forceSpline( f1, f1, r1, r2 ) {
  }

  // return ZBL energy (use with numerical diff)
  function energyZBL( Z1, Z2 ) {
    var a0 = 0.5,
        a = 0.8854 * a0 / ( Math.pow(Z1,0.23) + Math.pow(Z2,0.23) ),
        e0 = 1.0,
        pre = 1/(4*Math.PI*e0) * Z1*Z2 * e*e,
        A = [ 0.1818, 0.5099, 0.2802, 0.02817 ],
        B = [ -3.2, -0.9423, -0.4029, -0.2016 ];

    function calculateEnergy(r) {
      var phi = 0.0,
          i;
      for( i = 0; i < 4; ++i ) {
        phi += A[i]*Math.exp(B[i]*r/a);
      }
      return pre/r*phi;
    }

    return calculateEnergy;
  }

  //
  // export public interface
  //
  return {
    Atom : Atom,
    AtomType : AtomType,
    Barrier : Barrier,
    Simulation : Simulation,
    Vector : Vector,

    render : {
      atoms : renderAtoms,
      barriers : renderBarriers,
      forces : renderForces
    },
    force : {
      lennardJones : forceLJ,
      morse : forceMorse,
      tabulated : forceTabulated,
      diff : forceNumericalDiff
    }
  };
})();
