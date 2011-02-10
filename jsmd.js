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
    var dx = this.x - Math.round(this.x/w) * w;
    var dy = this.y - Math.round(this.y/h) * h;
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
    var dx = a.x-b.x; 
    dx -= Math.round(dx/w) * w;
    var dy = a.y-b.y;
    dy -= Math.round(dy/h) * h;
    return Math.sqrt( dx*dx + dy*dy );
  }
  Vector.pbcdistance2 = function(a,b,w,h) {
    var dx = a.x-b.x; 
    dx -= Math.round(dx/w) * w;
    var dy = a.y-b.y;
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
    var sv = new Vector( this.p[0].x - a.x, this.p[0].y - a.y );
    var s = -(sv.x*this.c.x + sv.y*this.c.y);
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
  function Simulation(w,h) {
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
    this.nl = { dr : 0.0, data : [] };

    // timestep data
    this.dt = 0.01;
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
    var rm = rc+rp;
    this.rm2 = rm*rm; // safe distance squared (for fast neighborlist builds)
    this.rp = rp;
    this.rc = rc;

    // number and size of cells
    this.lc.nx = Math.floor(this.w/rm);
    this.lc.ny = Math.floor(this.h/rm);
    this.lc.dx = this.w/this.lc.nx;
    this.lc.dy = this.h/this.lc.ny;

    // build the 2D linkcell grid
    var l = new Array(this.lc.nx);
    for( var i = 0; i < this.lc.nx; ++i ) {
      l[i] = new Array(this.lc.ny);
    }
    this.lc.data = l;
  }
  Simulation.prototype.clearLinkcell = function() {
    // clear each cell
    for( var i = 0; i < this.lc.nx; ++i ) {
      for( var j = 0; j < this.lc.ny; ++j ) {
        this.lc.data[i][j] = [];
      }
    }
  }
  Simulation.prototype.updateLinkcell = function() {
    // clear first
    this.clearLinkcell();

    // repopulate with all atoms
    var lx, ly;
    for( i = 0; i < this.atoms.length; ++i ) {
      lx = Math.floor(this.atoms[i].p.x/this.lc.dx) % this.lc.nx;
      ly = Math.floor(this.atoms[i].p.y/this.lc.dy) % this.lc.ny;
      this.lc.data[lx][ly].push(i);
    }
  }
  Simulation.prototype.clearNeighborlist = function() {
    for( i = 0; i < this.atoms.length; ++i ) {
      this.nl.data[i] = [];
    }
  }
  Simulation.prototype.updateNeighborlist = function(dr) {
    // check if update is necessary (call without parameter to force update)
    this.nl.dr += 2.0*dr;
    if( dr !== undefined && this.nl.dr < this.rp ) { 
      return;
    }

    // update Linkcells first
    this.updateLinkcell()

    // clear Neigborlist
    this.clearNeighborlist()

    var i,j,i2,j2,k,l,m;
    var ka, la;
    var n = [ [0,1],[1,this.h-1],[1,0],[1,1] ]; // half the neighbors
    // loop over all cells
    for( i = 0; i < this.lc.nx; ++i ) {
      for( j = 0; j < this.lc.ny; ++j ) {
        // loop over all atoms in local cell
        var ll = this.lc.data[i][j].length;
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
  }
  Simulation.prototype.updateForces = function() {
    var i,j,k;  // integer
    var F,dr; // float
    var f; // function
    var rvec; // Vector

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
        dr = rvec.len();
        if( dr < this.rc ) {
          f = this.interaction[this.atoms[i].t][this.atoms[j].t];
          if( f !== undefined ) {
            F = f( dr, [ this.atoms[i].t, this.atoms[j].t ] );
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
            F = f( dr, [ this.atoms[i].t, this.barriers[j].t ] );
            //$('#mdlog').text(dr+','+f);
            rvec.scale(F/dr);
            this.atoms[i].f.add(rvec);
          }
        }
      }
    }
  }
  Simulation.prototype.velocityVerlet = function() {
    var i,j;
    var dt = this.dt;
    var rmax = 0.0, vmax = 0.0, amax = 0.0;

    // first velocity verlet step
    var dp = new jsmd.Vector(), dp2;
    for( i = 0; i < this.atoms.length; ++i ) {
      var m = this.types[this.atoms[i].t].m;
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

      // maximum velocity and acceleration
      vmax = Math.max( vmax, this.atoms[i].v.len2() );
      amax = Math.max( amax, this.atoms[i].f.len2()/(0.25*m*m) );
    }

    // compute timestep
    var dmax = 0.01;
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
    var c = this.canvas.ctx;
    var i;

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

  // built-in render routines
  function renderAtoms(c) {
    function drawAtom(x,y,r) {
      c.beginPath();
      c.arc( x, y, r, 0, Math.PI*2.0, true);
      c.closePath();
      c.fill();
    }
    // draw nuclei
    c.strokeStyle = "rgba(0,100,0,0.5)";
    var x,y,r;
    for( var i = 0; i < this.atoms.length; ++i ) {
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
  function renderBarriers(c) {
    // draw barriers
    c.strokeStyle = "rgba(0,0,100,1)";
    for( var i = 0; i < this.barriers.length; ++i ) {
      c.beginPath();
      c.moveTo( this.barriers[i].p[0].x, this.barriers[i].p[0].y );
      c.lineTo( this.barriers[i].p[1].x, this.barriers[i].p[1].y );
      c.closePath();
      c.stroke();
    }
  }
  function renderForces(c) {
    // draw forces
    c.strokeStyle = "rgba(0,100,0,0.5)";
    c.lineWidth = 0.01;
    for( var i = 0; i < this.atoms.length; ++i ) {
      c.beginPath();
      c.moveTo( this.atoms[i].p.x, this.atoms[i].p.y );
      c.lineTo( this.atoms[i].p.x + this.atoms[i].f.x * 0.01, this.atoms[i].p.y + this.atoms[i].f.y * 0.01 );
      c.closePath();
      c.stroke();
    }
  }

  // built-in force functions
  function forceLJ(r,t) {
    r /= 4;
    var e = 10.0;
    var rm6 = 1.0;
    var rm12 = 1.0;
    return 12.0*e*( rm6*Math.pow(r,-7.0) - rm12*Math.pow(r,-13.0) );
  }
  function forceMorse(r,t) {
    //De*( 1-exp(-a*(x-re)) )**2, (x-re)**2
    var De = 1.0;
    var re = 1.0; // 0.5*(t[0].re+t[1].re)
    var a = 1.0; // Math.sqrt(ke/(2*De))
    var ex = Math.exp(-a*(r-re));
    return 2*a*De*(1-ex)*ex;
  }

  // export public interface
  return {
    Atom : Atom,
    AtomType : AtomType,
    Barrier : Barrier,
    Vector : Vector,
    Simulation : Simulation,

    render : {
      atoms : renderAtoms,
      barriers : renderBarriers,
      forces : renderForces
    },
    force : {
      lennardJones : forceLJ,
      morse : forceMorse
    }
  };
})();
