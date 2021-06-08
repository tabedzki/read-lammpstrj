#include "globals.h"
#include "timing.h"

void calc_Dweights( void ) ;
int remove_particles(int,int*,int,int*);
void swap_particles( void ) ;
void swap_ghosts( void ) ;
void particle_orientations( void ) ;
void particle_Stensor( void ) ;


void update_positions() {

  move_t_in = time(0) ;

  int i,j, id ;

  Nbound_ct = Sbound_ct = 0 ;

  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;

    double dsq = Diff[ mol_type[i] ] * delt ;


    // Move the particles, keep them in the box //
    for ( j=0 ; j<Dim ; j++ ) 
      x[id][j] = x[id][j] + dsq * f[id][j] + sqrt( 2.0 * dsq ) * gasdev2() ;


    if ( nprocs > 1 ) {
      // Track the particles that move to adjacent processors //
      if ( x[id][Dim-1] < z_min ) {
        send_S_inds[ Sbound_ct ] = id ;
        Sbound_ct++ ;
      }
      else if ( x[id][Dim-1] > z_max ) {
        send_N_inds[ Nbound_ct ] = id ;
        Nbound_ct++ ;
      }
    }

    for ( j=0 ; j<Dim ; j++ ) {
      if ( x[id][j] > L[j] ) 
        x[id][j] -= L[j] ;

      else if ( x[id][j] < 0.0 )
        x[id][j] += L[j] ;
    }

  }//for ( i=0:ns_loc
  

  if ( mu != 0.0 || print_vectors == 1 || print_order_params == 1) {
    particle_orientations() ;
  }



  ///////////////////////////////////////////////////////
  // Remove departing particles from current processor //
  ///////////////////////////////////////////////////////
  if ( nprocs > 1 ) {
    ns_loc = remove_particles( ns_loc, my_inds, Nbound_ct, send_N_inds ) ;
    ns_loc = remove_particles( ns_loc, my_inds, Sbound_ct, send_S_inds ) ;
  }

  move_t_out = time(0) ;
  move_minus_comm += (move_t_out - move_t_in) ;
 


  ////////////////////////////////////////////////////////
  // Send/receive the particles who are switching procs //
  ////////////////////////////////////////////////////////
  if ( nprocs > 1 ) {
    swap_particles() ;

    swap_ghosts() ;
  }


  // Calculate particle S tensors
  if ( mu != 0.0 ) {
    particle_Stensor() ;
  }



  move_t_out = time(0) ;
  move_tot_time += (move_t_out - move_t_in) ;

}


