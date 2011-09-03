function Neighborlist2d(lc) {
  this.lc = lc;
  this.dr = 0.0;
  this.data = null;
  this.count = 0;
}

Neighborlist2d.prototype.clear = function() {
  this.data = [];
  for( var i = 0; i < this.lc.sim.atoms.length; ++i ) {
    this.data[i] = [];
  }
}
Neighborlist2d.prototype.update = function(dr) {
  var sim = this.lc.sim;
  // check if update is necessary (call without parameter to force update)
  this.dr += 2.0*dr;
  if( this.dr < sim.rp && 
      this.data !== null &&
      this.data.length == this.lc.sim.atoms.length ) { 
    return;
  }
  this.dr = 0.0;

  // update Linkcells first
  this.lc.update()

  // clear Neigborlist
  this.clear()

  var i,j,i2,j2,n,l,m,
      ka, la, ll,
      // half the neighbors
      neigh = [ [0,1],[1,this.lc.ny-1],[1,0],[1,1] ];

  // loop over all cells
  for( i = 0; i < this.lc.nx; ++i ) {
    for( j = 0; j < this.lc.ny; ++j ) {

      // loop over all atoms in local cell
      ll = this.lc.data[i][j].length;
      for( n = 0; n < ll; ++n ) {
        // loop over remaining atoms in local cell
        ka = this.lc.data[i][j][n];
        for( l = n+1; l < ll; ++l ) {
          la = this.lc.data[i][j][l];
          if( Vector2d.pbcdistance2( sim.atoms[ka].p, sim.atoms[la].p, sim.ss ) < this.lc.rm2 ) {
            this.data[ka].push(la);
          }
        }

        // loop over half the neighbor cells
        for( m = 0; m < neigh.length; ++m ) {
          i2 = (i+neigh[m][0]) % this.lc.nx;
          j2 = (j+neigh[m][1]) % this.lc.ny;
          // loop over all atoms in those neighbor cells
          for( l = 0; l < this.lc.data[i2][j2].length; ++l ) {
            la = this.lc.data[i2][j2][l];
            if( Vector2d.pbcdistance2( sim.atoms[ka].p, sim.atoms[la].p, sim.ss ) < this.lc.rm2 ) {
              this.data[ka].push(la);
            }
          }
        }
      }
      // end local cell

    }
  }

  // increase update counter
  this.count++;
}
