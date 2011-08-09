function initJSMD(dim) {
  // depending on which vector class was loaded
  var Vector, Linkcell, Neighborlist;
            
  // Atom constructor
  function Atom(x,y,z) {
    this.p = new Vector(x,y,z); // y may be 'undefined' if x is a vector!
    this.v = new Vector();
    this.f = new Vector();
    this.t = 0;
  }

  // AtomType constructor
  function AtomType() {
    this.r = 3;
    this.color  = "rgb(255,0,0)";
    this.Z = 1;
    this.m = 1.0;
  }

  // Barrier constructor strictly 2d right now
  function Barrier(p1,p2,sim) {
    this.p = [ new Vector(p1),  new Vector(p2) ];

    // normal vector
    this.n = Vector.sub(p2,p1);
    this.n.normalize();

    // colinear vector
    this.c = new Vector( p2.x-p1.x, p2.y-p1.y );
    this.l = this.c.len(); // length
    this.c.scale(1.0/this.l);

    // store force acting on the barrier
    this.f = new Vector( 0.0, 0.0 );

    // default interation type
    this.t = 0;
    
    // link back to ssimulation object for PBC
    this.sim = sim;
  }
  Barrier.prototype.dist = function(a) {
    var sv = new Vector( this.p[0].x - a.x, this.p[0].y - a.y ), s;

    sv.dwrap( this.sim.ss );
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

  // helper function to merge arguments with defaults
  function mergeArgs( arg, def ) {
    // no arguments supplied, just return defaults
    if( arg === undefined ) {
      return def;
    }
    // set unset options in arg to the defaults
    for( i in def ) {
      if( def.hasOwnProperty(i) ) {}
    }
  }

  // Simulation constructor (this contains all the important logic)
  function Simulation(ss, options) {
    this.atoms    = []; // list of atoms
    this.barriers = []; // list of barriers
    this.types    = []; // list of atom types

    // interaction matrix
    this.interaction = [];

    // box dimensions
    this.ss = new Vector();
    this.ss.set(ss);

    // canvas for visualization
    this.canvas = {};

    // setup default render chain
    this.renderChain = [ renderAtoms, renderForces, renderBarriers ];

    // setup default compute chain
    this.computeChain = [ computeVerlet1, computeWrapUpdate, computeForces, computeVerlet2 ];
    
    // timestep data
    this.dt = 0.0;   // set by dynamic timestepper
    this.step = 0;   // number of current MD step
    this.time = 0.0; // expired simulation time

    // maximum spatial step for dynamic timestep algorithm
    this.dmax = 0.025;

    // drag factor (set to 1 to disable drag)
    this.drag = 0.995;
    
    // gravitation/wind
    this.f = new Vector(0,0);
    
    // store the virial here (for pressure calculation)
    this.vir = 0.0;
    // instantaneous temperature 
    this.T = 0.0;
    // instantaneous pressure
    this.P = 0.0;
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

    // initialize Linkcell and Neighborlist classes
    this.lc = new Linkcell(this);
    this.nl = new Neighborlist(this.lc);
  }
 
  function computeForces(store) {
    var i,j,k,  // integer
        F,dr, // float
        f, // function
        rvec; // Vector

    // zero forces set timestep
    for( i = 0; i < this.atoms.length; ++i ) {
      this.atoms[i].f.set(this.f);
    }
    for( i = 0; i < this.barriers.length; ++i ) {
      this.barriers[i].f.zero();
    }
    
    // clear virial
    this.vir = 0.0;

    // build up forces (newtonian)
    for( i = 0; i < this.atoms.length; ++i ) {
      //for( j = i+1; j < atoms.length; ++j ) {     // old iteration over all other atoms
      for( k = 0; k < this.nl.data[i].length; ++k ) {  // new iteration over neighborlist
        j = this.nl.data[i][k];

        rvec = jsmd.Vector.sub( this.atoms[j].p, this.atoms[i].p);
        rvec.dwrap(this.ss);
        dr = rvec.len2();
        if( dr < this.rc2 ) {
          dr = Math.sqrt(dr);
          f = this.interaction[this.atoms[i].t][this.atoms[j].t];
          if( f !== undefined ) {
            F = f.call( this, dr, [ this.atoms[i].t, this.atoms[j].t ] );
            rvec.scale(F/dr);
            this.atoms[i].f.add(rvec);
            this.atoms[j].f.sub(rvec);
            
            this.vir += F*dr;
          }
        }
      }
    }

    // add barrier interaction
    for( i = 0; i < this.atoms.length; ++i ) {
      for( j = 0; j < this.barriers.length; ++j ) {
        // find distance to barrier
        rvec = this.barriers[j].dist(this.atoms[i].p);
        dr = rvec.pbclen(this.ss);
        if( dr < this.rc ) {
          f = this.interaction[this.atoms[i].t][this.barriers[j].t];
          if( f !== undefined ) {
            F = f.call( this, dr, [ this.atoms[i].t, this.barriers[j].t ] );
            console.log(dr+','+f);
            rvec.scale(F/dr);
            this.atoms[i].f.add(rvec);
            this.barriers[j].f.sub(rvec);
            
            this.vir += F*dr;
          }
        }
      }
    }
    
    // virial (1/3.0 in 3d, 1/2.0 in 2d)
    this.vir /= 2.0;
  }
  
  // first velocity verlet step
  function computeVerlet1(store) {
    var i, m, v2, 
        dp = new jsmd.Vector();
    store.rmax = 0.0;

    // first velocity verlet step
    for( i = 0; i < this.atoms.length; ++i ) {
      m = this.types[this.atoms[i].t].m;
      dp.set( this.atoms[i].v ); dp.scale(this.dt); // dp = v*dt
      dp.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*this.dt*this.dt) );
      this.atoms[i].p.add(dp);
      this.atoms[i].v.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*this.dt) );

      // track maximum displacement
      store.rmax = Math.max( store.rmax, dp.len2() );
    }
  }
  
  // wrap atoms and update neighbor lists (needs store.rmax)
  function computeWrapUpdate(store) {
    var i;
    for( i = 0; i < this.atoms.length; ++i ) {
      this.atoms[i].p.wrap(this.ss);
    }
    this.nl.update(Math.sqrt(store.rmax));
  }
  
  // wrap atoms and update neighbor lists (needs store.rmax)
  function computeBounce(store) {
    var i, p, v;
    for( i = 0; i < this.atoms.length; ++i ) {
      p = this.atoms[i].p;
      v = this.atoms[i].v;
      if( p.x < 0 ) {
        p.x = -p.x; v.x = -v.x;
      } else if ( p.x > this.ss.x ) {
        p.x = 2.0*this.ss.x - p.x; v.x = -v.x;
      }
      if( p.y < 0 ) {
        p.y = -p.y; v.y = -v.y;
      } else if ( p.y > this.ss.y ) {
        p.y = 2.0*this.ss.y - p.y; v.y = -v.y;
      }
      if( p.z < 0 ) {
        p.z = -p.z; v.z = -v.z;
      } else if ( p.z > this.ss.z ) {
        p.z = 2.0*this.ss.z - p.z; v.z = -v.z;
      }
    }
    //this.nl.update(Math.sqrt(store.rmax));
  }
  
  // second velocity verlet step
  function computeVerlet2(store) {
    var i, m, v2, vmax = 0.0, amax = 0.0, dmax;
    
    this.T = 0.0;
    for( i = 0; i < this.atoms.length; ++i ) {
      m = this.types[this.atoms[i].t].m;
      this.atoms[i].v.add( jsmd.Vector.scale(this.atoms[i].f, 0.5/m*this.dt) );

      // linear drag term
      this.atoms[i].v.scale(this.drag);

      // calculate temperature
      v2 = this.atoms[i].v.len2();
      this.T += this.types[this.atoms[i].t].m * v2;
      
      // maximum velocity and acceleration
      vmax = Math.max( vmax, v2 );
      amax = Math.max( amax, this.atoms[i].f.len2()/(0.25*m*m) );
    }
    this.T /= 2.0; // Vector.dim
    
    // calculate pressure (PV=NkBT-this.vir)
    this.P = ( this.T -this.vir ) / ( this.ss.vol() );
    
    // 3/2*N*kB*T = 1/2*sum(m*v^2) (2/2NkT in 2d?)
    this.T = this.T/(this.atoms.length); // *1/kB

    // increase step counters
    this.step++;
    this.time += this.dt;

    // compute timestep
    vmax = Math.sqrt(vmax);
    amax = Math.sqrt(amax);
    this.dt = Math.min( 0.01, this.dmax/vmax, Math.sqrt(2*this.dmax/amax) );
  }

  // Berendsen hydrostatic barostat factory function
  function computeBerendsenP( P0, tau ) {
    return function(store) {
      var i, l = ( 1.0 - this.dt/tau * ( P0 - this.P ) );

      // scale atomic coordinates
      for( i = 0; i < this.atoms.length; ++i ) {
        this.atoms[i].p.scale(l);
      }
      
      // scale box
      this.ss.scale(l);
    }
  }
  
  // run one full timestep, process all items in the compute chain
  Simulation.prototype.run = function(steps) {
    // hash object passed as reference to allow datatransfer between computeChain members
    var store = {}, i, j;

    // run multiple steps
    for( i = 0; i < steps; ++i ) {
      // process compute chain
      for( j = 0; j < this.computeChain.length; ++j ) {
        this.computeChain[j].call(this,store);
      }
    }
    
    return store;
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
    c.setTransform( this.canvas.w/this.ss.x, 0,0, this.canvas.h/this.ss.y, 0,0 );

    // clear
    c.fillStyle = "rgb(200,200,200)";
    c.fillRect(0, 0, this.ss.x,this.ss.y );

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
      if( x <= r || y <= r || x+r > this.ss.x || y+r > this.ss.y ) {
        if( x <= r ) drawAtom(x+this.ss.x,y,r);
        if( y <= r ) drawAtom(x,y+this.ss.y,r);
        if( x <= r &&  y <= r ) drawAtom(x+this.ss.x,y+this.ss.y,r);
        if( x+r > this.ss.x ) drawAtom(x-this.ss.x,y,r);
        if( y+r > this.ss.y ) drawAtom(x,y-this.ss.y,r);
        if( x+r > this.ss.x && y+r > this.ss.y) drawAtom(x-this.ss.x,y-this.ss.y,r);
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

    return function(r) {
      r /= 4;
      return 12.0*e*( rm6*Math.pow(r,-7.0) - rm12*Math.pow(r,-13.0) );
    }
  }

  // hard wall
  function force12() {
    return function(r) {
      return -Math.pow(r,-13.0);
    }
  }

  // Morse potential for pair equilibrium distance re, potential softness 1/a, and well depth De
  function forceMorse( re_, a_, De_ ) {
    //De*( 1-exp(-a*(x-re)) )**2, (x-re)**2
    var re = re_ !== undefined ? re_ : 1.5,
        a  = a_  !== undefined ? a_  : 2.0,
        De = De_ !== undefined ? De_ : 1.0;

    return function(r) {
      var ex = Math.exp( -a*(r-re) );
      return 2.0 * a * De * (1-ex) * ex;
    }
  }

  // Morse potential energy
  function energyMorse( re_, a_, De_ ) {
    //De*( 1-exp(-a*(x-re)) )**2, (x-re)**2
    var re = re_ !== undefined ? re_ : 1.5,
        a  = a_  !== undefined ? a_  : 2.0,
        De = De_ !== undefined ? De_ : 1.0;

    return function(r) {
      var ex = Math.exp( -a*(r-re) );
      return De * (1-ex) * (1-ex);
    }
  }

  // return a force linearly interpolated from tabulated values of function f (fast!)
  function forceTabulated(f, dr_, rc_, noint_ ) {
    // parameters and pretabulated values stored in closure
    var dr = dr_ !== undefined ? dr_ : 0.01,
        rc = rc_ !== undefined ? rc_ : 10.0,
        doint = noint !== true,
        mul = rc/dr,
        table, n = Math.round(mul*rc+20),
        i;
        
    if (typeof Float32Array != "undefined") {   
      table = new Float32Array(n);
    } else {
      table = new Array(n);
    }

    for( i = 0; i <= mul*rc+10; ++i ) {
      table[i] = f.call( this, i/mul );
    }

    // linear interpolation function
    function linint(r) {
      r *= mul;
      var b = Math.floor(r);
      r -= b;
      return (1-r)*table[b] + r*table[b+1];
    }
    // no interpolation
    function noint(r) {
      return table[Math.floor(r*mul)];
    }
    
    // return evaluation function (called from updateForces)
    return doint ? linint : noint;
  }

  // return numerical derivative of f
  function forceNumericalDiff( f, dr_ ) {
    var dr = dr_ !== undefined ? dr_ : 0.0001;
    
    return function (r) {
       return 0.5 * ( f.call(this,r+dr) - f.call(this,r-dr) ) / dr;
    }
  }

  // ZBL fore function for nuclear charges Z1 and Z2
  function forceZBL( Z1_, Z2_ ) {
  }

  // return force splined together from f1 for r<=r1, a spline for r1<r<=r2, and f2 for r2<r
  function forceSpline( f1, f2, r1, r2 ) {
  }

  // return ZBL energy (use with numerical diff)
  function energyZBL( Z1, Z2 ) {
    var a0 = 0.539177,
        a = 0.8854 * a0 / ( Math.pow(Z1,0.23) + Math.pow(Z2,0.23) ),
        e = 1.0,
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

  // configure dimension
  var missing = [];
  if( dim === 2 ) {
    if( window.Vector2d === undefined ) { missing.push('vector2d.js'); }
    if( window.Linkcell2d === undefined ) { missing.push('linkcell2d.js'); }
    if( window.Neighborlist2d === undefined ) { missing.push('neighborlist2d.js'); }
    if( missing.length > 0 ) {
      alert('Must load ' + missing.join(',') );
    }
    Vector = window.Vector2d;
    Linkcell = window.Linkcell2d;
    Neighborlist = window.Neighborlist2d;
  } else if( dim === 3 ) {
    if( window.Vector3d === undefined ) { missing.push('vector3d.js'); }
    if( window.Linkcell3d === undefined ) { missing.push('linkcell3d.js'); }
    if( window.Neighborlist3d === undefined ) { missing.push('neighborlist3d.js'); }
    if( missing.length > 0 ) {
      alert('Must load ' + missing.join(',') );
    }
    Vector = window.Vector3d;
    Linkcell = window.Linkcell3d;
    Neighborlist = window.Neighborlist3d;
  }
  else {
    alert('Dimension parameter must be 2 or 3');
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

    compute : {
      forces : computeForces,
      verlet1 : computeVerlet1,
      verlet2 : computeVerlet2,
      bounce : computeBounce,
      wrapUpdate : computeWrapUpdate,
      berendsenP: computeBerendsenP
    },
    render : {
      atoms : renderAtoms,
      barriers : renderBarriers,
      forces : renderForces
    },
    force : {
      lennardJones : forceLJ,
      wall : force12,
      morse : forceMorse,
      tabulated : forceTabulated,
      diff : forceNumericalDiff
    },
    energy : {
      morse : energyMorse,
      ZBL : energyZBL
    }
  };
};
