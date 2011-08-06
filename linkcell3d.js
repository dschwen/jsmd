function Linkcell3d( sim ) {
  var l,i,j,rm;
  
  // link to host simulation object
  this.sim = sim;

  // cut off must be set
  rm = sim.rc + sim.rp;
  this.rm2 = rm*rm;

  // number and size of cells (have at least 3 cells in each direction)
  this.nx = Math.max( 3, Math.floor(sim.ss.x/rm) );
  this.ny = Math.max( 3, Math.floor(sim.ss.y/rm) );
  this.nz = Math.max( 3, Math.floor(sim.ss.z/rm) );
  this.dx = sim.ss.x/this.nx;
  this.dy = sim.ss.y/this.ny;
  this.dz = sim.ss.z/this.nz;

  // build the 3D linkcell grid
  l = new Array(this.nx);
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

  // repopulate with all atoms
  var lx, ly,lz, i;
  for( i = 0; i < this.sim.atoms.length; ++i ) {
    lx = Math.floor(this.sim.atoms[i].p.x/this.dx) % this.nx;
    ly = Math.floor(this.sim.atoms[i].p.y/this.dy) % this.ny;
    lz = Math.floor(this.sim.atoms[i].p.z/this.dz) % this.nz;
    this.data[lx][ly][lz].push(i);
  }
}
