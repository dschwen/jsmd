function Linkcell2d( sim ) {
  var l,i,rm;
  
  // link to host simulation object
  this.sim = sim;

  // cut off must be set
  rm = sim.rc + sim.rp;
  this.rm2 = rm*rm;

  // number and size of cells (have at least 3 cells in each direction)
  this.nx = Math.max( 3, Math.floor(sim.ss.x/rm) );
  this.ny = Math.max( 3, Math.floor(sim.ss.y/rm) );
  this.dx = sim.ss.x/this.nx;
  this.dy = sim.ss.y/this.ny;

  // build the 2D linkcell grid
  l = new Array(this.nx);
  for( i = 0; i < this.nx; ++i ) {
    l[i] = new Array(this.ny);
  }
  this.data = l;
}

Linkcell2d.prototype.clear = function() {
  // clear each cell
  var i,j;
  for( i = 0; i < this.nx; ++i ) {
    for( j = 0; j < this.ny; ++j ) {
      this.data[i][j] = [];
    }
  }
}

Linkcell2d.prototype.update = function() {
  // clear first
  this.clear();

  // repopulate with all atoms
  var lx, ly, i;
  for( i = 0; i < this.sim.atoms.length; ++i ) {
    lx = Math.floor(this.sim.atoms[i].p.x/this.dx) % this.nx;
    ly = Math.floor(this.sim.atoms[i].p.y/this.dy) % this.ny;
    this.data[lx][ly].push(i);
  }
}
