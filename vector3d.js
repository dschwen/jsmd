// Vector constructor 
function Vector3d(x,y,z) { // or Vector2d(v), where v is a Vector, or Vector2d(), which initializes to zero
  if( y === undefined && x !== undefined ) {
    this.x = x.x;
    this.y = x.y;
    this.z = x.z;
  } else {
    this.x = x || 0.0;
    this.y = y || 0.0;
    this.z = z || 0.0;
  }
}
Vector3d.prototype.len = function() {
  return Math.sqrt( this.x * this.x + this.y * this.y + this.z * this.z );
}
Vector3d.prototype.len2 = function() {
  return this.x * this.x + this.y * this.y + this.z * this.z;
}
Vector3d.prototype.pbclen = function(ss) {
  var dx = this.x - Math.round(this.x/ss.x) * ss.x,
      dy = this.y - Math.round(this.y/ss.y) * ss.y,
      dz = this.z - Math.round(this.z/ss.z) * ss.z;
  return Math.sqrt( dx*dx + dy*dy + dz*dz );
}
Vector3d.prototype.normalize = function() {
  var len = this.len();
  if( len === 0 ) {
    this.x = 1; this.y = 0; this.z = 0;
  } else {
    this.x /= len; this.y /= len; this.z /= len;
  }
}
Vector3d.prototype.vol = function() {
  return this.x * this.y * this.z;
}
Vector3d.distance = function(a,b) {
  return Math.sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
}
Vector3d.distance2 = function(a,b) {
  return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
}
Vector3d.pbcdistance = function(a,b,ss) {
  var dx = a.x-b.x,
      dy = a.y-b.y,
      dz = a.z-b.z;
  dx -= Math.round(dx/ss.x) * ss.x;
  dy -= Math.round(dy/ss.y) * ss.y;
  dz -= Math.round(dz/ss.z) * ss.z;
  return Math.sqrt( dx*dx + dy*dy + dz*dz );
}
Vector3d.pbcdistance2 = function(a,b,ss) {
  var dx = a.x-b.x,
      dy = a.y-b.y,
      dz = a.z-b.z;
  dx -= Math.round(dx/ss.x) * ss.x;
  dy -= Math.round(dy/ss.y) * ss.y;
  dz -= Math.round(dz/ss.z) * ss.z;
  return dx*dx + dy*dy + dz*dz;
}
Vector3d.prototype.wrap = function(ss) {
  this.x -= Math.floor(this.x/ss.x) * ss.x;
  this.y -= Math.floor(this.y/ss.y) * ss.y;
  this.z -= Math.floor(this.z/ss.z) * ss.z;
}
Vector3d.prototype.dwrap = function(ss) {
  this.x -= Math.round(this.x/ss.x) * ss.x;
  this.y -= Math.round(this.y/ss.y) * ss.y;
  this.z -= Math.round(this.z/ss.z) * ss.z;
}
Vector3d.sub = function(a,b) {
  return new Vector3d( a.x-b.x, a.y-b.y, a.z-b.z );
}
Vector3d.prototype.sub = function(a) {
  this.x -= a.x;
  this.y -= a.y;
  this.z -= a.z;
}
Vector3d.add = function(a,b) {
  return new Vector3d( a.x+b.x, a.y+b.y, a.z+b.z );
}
Vector3d.prototype.add = function(a) {
  this.x += a.x;
  this.y += a.y;
  this.z += a.z;
}
Vector3d.scale = function(a,b) {
  return new Vector3d( a.x*b, a.y*b, a.z*b );
}
Vector3d.prototype.scale = function(a) {
  this.x *= a;
  this.y *= a;
  this.z *= a;
}
Vector3d.smul = function(a,b) { // scalar product
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vector3d.prototype.smul = function(a) { // scalar product
  return this.x*a.x + this.y*a.y + this.z*a.z;
}
Vector3d.prototype.proj = function(a) { // assume a is a unit vector!
  var s = this.smul(a);
  this.x = a.x * s;
  this.y = a.y * s;
  this.z = a.z * s;
}
Vector3d.prototype.zero = function() {
  this.x = 0.0;
  this.y = 0.0;
  this.z = 0.0;
}
Vector3d.prototype.set = function(a) {
  this.x = a.x;
  this.y = a.y;
  this.z = a.z;
}
Vector3d.random = function(ss) {
  return new Vector3d( Math.random()*ss.x, Math.random()*ss.y, Math.random()*ss.z );
}

