function Linkcell2d( sim ) {
  // link to host simulation object
  this.sim = sim;

  // cut off must be set
  var rm = sim.rc + sim.rp;
  this.rm2 = rm*rm;

  // number and size of cells (have at least 3 cells in each direction)
  this.nx = 3;
  this.ny = 3;

  // build the 2D linkcell grid
  this.setup();
}

Linkcell2d.prototype.setup = function() {
  // build the 2D linkcell grid
  var i, j, l = new Array(this.nx);
  for( i = 0; i < this.nx; ++i ) {
    l[i] = new Array(this.ny);
    for( j = 0; j < this.ny; ++j ) {
      l[i][j] = [];
    }
  }
  this.data = l;
}

Linkcell2d.prototype.clear = function() {
  // clear each cell
  var i,j;
  for( i = 0; i < this.nx; ++i ) {
    for( j = 0; j < this.ny; ++j ) {
      //this.data[i][j] = [];
      this.data[i][j].length = 0;
    }
  }
}

Linkcell2d.prototype.update = function() {

  // calculate optimal linkcell size (box size can be dynamic)
  var rm = sim.rc + sim.rp,
      nnx = Math.max( 3, Math.floor(sim.ss.x/rm) ),
      nny = Math.max( 3, Math.floor(sim.ss.y/rm) );

  // do we need to redimension the linkcell grid?
  if( this.nx != nnx || this.ny != nny ) {
    this.nx = nnx;
    this.ny = nny;
    this.setup();
  } else {
    // clear
    this.clear();
  }

  // cell size
  var dx = sim.ss.x/this.nx,
      dy = sim.ss.y/this.ny;

  // repopulate with all atoms
  var lx, ly, i;
  for( i = 0; i < this.sim.atoms.length; ++i ) {
    lx = Math.floor(this.sim.atoms[i].p.x/dx) % this.nx;
    ly = Math.floor(this.sim.atoms[i].p.y/dy) % this.ny;
    this.data[lx][ly].push(i);
  }
}
