function Linkcell3d( sim ) {
  var j,rm;
  
  // link to host simulation object
  this.sim = sim;

  // cut off must be set
  rm = sim.rc + sim.rp;
  this.rm2 = rm*rm;

  // number and size of cells (have at least 3 cells in each direction)
  this.nx = 3;
  this.ny = 3;
  this.nz = 3;

  // build the 3D linkcell grid
  this.setup();
}

Linkcell3d.prototype.setup = function() {
  // build the 3D linkcell grid
  var i, l = new Array(this.nx);
  for( i = 0; i < this.nx; ++i ) {
    l[i] = new Array(this.ny);
    for( j = 0; j < this.ny; ++j ) {
      l[i][j] = new Array(this.nz);
    }
  }
  this.data = l;
}

Linkcell3d.prototype.clear = function() {
  // clear each cell
  var i,j,k;
  for( i = 0; i < this.nx; ++i ) {
    for( j = 0; j < this.ny; ++j ) {
      for( k = 0; k < this.nz; ++k ) {
        this.data[i][j][k] = [];
      }
    }
  }
}

Linkcell3d.prototype.update = function() {
  // clear first
  this.clear();

  // calculate optimal linkcell size (box size can be dynamic)
  var rm = sim.rc + sim.rp,
      nnx = Math.max( 3, Math.floor(sim.ss.x/rm) ),
      nny = Math.max( 3, Math.floor(sim.ss.y/rm) ),
      nnz = Math.max( 3, Math.floor(sim.ss.z/rm) );

  // do we need to redimension the linkcell grid?
  if( this.nx != nnx || this.ny != nny || this.nz != nnz ) {
    this.nx = nnx;
    this.ny = nny;
    this.nz = nnz;
    this.setup();
  }
   
  // cell size
  var dx = sim.ss.x/this.nx,
      dy = sim.ss.y/this.ny,
      dz = sim.ss.z/this.nz;

  // repopulate with all atoms
  var lx, ly,lz, i;
  for( i = 0; i < this.sim.atoms.length; ++i ) {
    lx = Math.floor(this.sim.atoms[i].p.x/dx) % this.nx;
    ly = Math.floor(this.sim.atoms[i].p.y/dy) % this.ny;
    lz = Math.floor(this.sim.atoms[i].p.z/dz) % this.nz;
    this.data[lx][ly][lz].push(i);
  }
  
  // TODO: density (maximum cell occupation) could also be used as a hint for redimensioning the grid (gas!)
}
