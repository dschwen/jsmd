// Vector constructor 
function Vector2d(x,y) { // or Vector2d(v), where v is a Vector, or Vector2d(), which initializes to zero
  if( y === undefined && x !== undefined ) {
    this.x = x.x;
    this.y = x.y;
  } else {
    this.x = x || 0.0;
    this.y = y || 0.0;
  }
}
Vector2d.prototype.len = function() {
  return Math.sqrt( this.x * this.x + this.y * this.y );
}
Vector2d.prototype.len2 = function() {
  return this.x * this.x + this.y * this.y;
}
Vector2d.prototype.pbclen = function(ss) {
  var dx = this.x - Math.round(this.x/ss.x) * ss.x,
      dy = this.y - Math.round(this.y/ss.y) * ss.y;
  return Math.sqrt( dx*dx + dy*dy );
}
Vector2d.prototype.normalize = function() {
  var len = this.len();
  if( len === 0 ) {
    this.x = 1; this.y = 0;
  } else {
    this.x /= len; this.y /= len;
  }
}
Vector2d.prototype.vol = function() {
  return this.x * this.y;
}
Vector2d.distance = function(a,b) {
  return Math.sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) );
}
Vector2d.distance2 = function(a,b) {
  return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}
Vector2d.pbcdistance = function(a,b,ss) {
  var dx = a.x-b.x,
      dy = a.y-b.y;
  dx -= Math.round(dx/ss.x) * ss.x;
  dy -= Math.round(dy/ss.y) * ss.y;
  return Math.sqrt( dx*dx + dy*dy );
}
Vector2d.pbcdistance2 = function(a,b,ss) {
  var dx = a.x-b.x,
      dy = a.y-b.y;
  dx -= Math.round(dx/ss.x) * ss.x;
  dy -= Math.round(dy/ss.y) * ss.y;
  return dx*dx + dy*dy;
}
Vector2d.prototype.wrap = function(ss) {
  this.x -= Math.floor(this.x/ss.x) * ss.x;
  this.y -= Math.floor(this.y/ss.y) * ss.y;
}
Vector2d.prototype.dwrap = function(ss) {
  this.x -= Math.round(this.x/ss.x) * ss.x;
  this.y -= Math.round(this.y/ss.y) * ss.y;
}
Vector2d.sub = function(a,b) {
  return new Vector2d( a.x-b.x, a.y-b.y );
}
Vector2d.prototype.sub = function(a) {
  this.x -= a.x;
  this.y -= a.y;
}
Vector2d.add = function(a,b) {
  return new Vector2d( a.x+b.x, a.y+b.y );
}
Vector2d.prototype.add = function(a) {
  this.x += a.x;
  this.y += a.y;
}
Vector2d.scale = function(a,b) {
  return new Vector2d( a.x*b, a.y*b );
}
Vector2d.prototype.scale = function(a) {
  this.x *= a;
  this.y *= a;
}
Vector2d.smul = function(a,b) { // scalar product
  return a.x*b.x + a.y*b.y;
}
Vector2d.prototype.smul = function(a) { // scalar product
  return this.x*a.x + this.y*a.y;
}
Vector2d.prototype.proj = function(a) { // assume a is a unit vector!
  var s = this.smul(a);
  this.x = a.x * s;
  this.y = a.y * s;
}
Vector2d.prototype.zero = function() {
  this.x = 0.0;
  this.y = 0.0;
}
Vector2d.prototype.set = function(a) {
  this.x = a.x;
  this.y = a.y;
}
Vector2d.random = function(ss) {
  return new Vector2d( Math.random()*ss.x, Math.random()*ss.y );
}
Vector2d.prototype.random = function(ss) {
  this.x = Math.random()*ss.x;
  this.y = Math.random()*ss.y;
}

