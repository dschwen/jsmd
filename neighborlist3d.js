function Neighborlist3d(lc) {
  this.lc = lc;
  this.dr = 0.0;
  this.data = null;
  this.count = 0;
  this.pbc = true;
}

Neighborlist3d.prototype.clear = function() {
  this.data = [];
  for( var i = 0; i < this.lc.sim.atoms.length; ++i ) {
    this.data[i] = [];
  }
}
Neighborlist3d.prototype.update = function(dr) {
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

  var i,j,k,i2,j2,k2,n,l,m,
      ka, la, ll,
      // half the neighbors
      neigh = [ 
        [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1], 
        [this.lc.nx-1,1,0], [1,0,this.lc.nz-1], [0,this.lc.ny-1,1], 
        [this.lc.nx-1,1,1], [1,this.lc.ny-1,1],[ 1,1,this.lc.nz-1] 
      ];

  // loop over all cells
  for( i = 0; i < this.lc.nx; ++i ) {
    for( j = 0; j < this.lc.ny; ++j ) {
      for( k = 0; k < this.lc.nz; ++k ) {

        // loop over all atoms in local cell
        ll = this.lc.data[i][j][k].length;
        for( n = 0; n < ll; ++n ) {
          // loop over remaining atoms in local cell
          ka = this.lc.data[i][j][k][n];
          for( l = n+1; l < ll; ++l ) {
            la = this.lc.data[i][j][k][l];
            if( Vector3d.pbcdistance2( sim.atoms[ka].p, sim.atoms[la].p, sim.ss ) < this.lc.rm2 ) {
              this.data[ka].push(la);
            }
          }

          // loop over half the neighbor cells
          for( m = 0; m < neigh.length; ++m ) {
            i2 = (i+neigh[m][0]);
            if( this.pbc ) {
              i2 %= this.lc.nx;
            } else if( i2 >= this.lc.nx ) continue; //TODO: wrong -1 at 0 -> nx-1!!!!
	
            j2 = (j+neigh[m][1]);
            if( this.pbc ) {
              j2 %= this.lc.ny;
            } else if( j2 >= this.lc.ny ) continue;

            k2 = (k+neigh[m][2]);
            if( this.pbc ) {
              k2 %= this.lc.nz;
            } else if( k2 >= this.lc.nz ) continue;

            // loop over all atoms in those neighbor cells
            for( l = 0; l < this.lc.data[i2][j2][k2].length; ++l ) {
              la = this.lc.data[i2][j2][k2][l];
              if( Vector3d.pbcdistance2( sim.atoms[ka].p, sim.atoms[la].p, sim.ss ) < this.lc.rm2 ) {
                this.data[ka].push(la);
              }
            }
          }
        }
        // end local cell

      }
    }
  }

  // increase update counter
  this.count++;
}
