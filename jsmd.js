var jsmd = (function(){
  // Vector constructor
  function Vector(x,y) {
    this.x = x || 0.0;
    this.y = y || 0.0;
  }
  Vector.prototype.len = function() {
    return Math.sqrt( this.x * this.x + this.y * this.y );
  }
  Vector.prototype.normalize = function() {
    var len = this.len();
    if( len === 0 ) {
      this.x = 1; this.y = 0;
    } else {
      this.x /= len; this.y /= len;
    }
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
    this.l = this.c.len(); // length
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

  // export public interface
  return {
    Atom : Atom,
    AtomType : AtomType,
    Barrier : Barrier,
    Vector : Vector
  };
})();