#include "globals.h"
#include "timing.h"

void prepare_ghost_index_lists( void ) ;
void prepare_ghost_force_lists(void);
int mpi_Recv_2d( int, int, double**, int, int, MPI_Comm ) ;
int mpi_Send_2d( int, int, double**, int, int, MPI_Comm ) ;



void swap_particles() {

  swap_t_in = time(0);

  int i, send_proc, send_to, rec_from ;

  if ( nprocs == 1 )
    return ;

  // Send the particles south //
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc - 1 ;
    if ( send_to < 0 ) send_to = nprocs - 1 ;

    if ( send_proc == myrank ) {
      extract_from_array( Sbound_ct, send_S_inds, send_S_x, x ) ;
      for ( i=0 ; i<Sbound_ct ; i++ ) local_flag[ send_S_inds[i] ] = 0 ;

      MPI_Send( &Sbound_ct, 1, MPI_INT, send_to, 606+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_S_inds, Sbound_ct, MPI_INT, send_to, 707+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Sbound_ct, Dim, send_S_x, send_to, 808+send_to, MPI_COMM_WORLD ) ;

      // Send orientation vectors if necessary
      if ( mu != 0.0  ) {
        extract_from_array( Sbound_ct, send_S_inds, send_S_x, mono_u ) ;
        mpi_Send_2d( Sbound_ct, Dim, send_S_x, send_to, 828+send_to, MPI_COMM_WORLD ) ;
      }
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_S_ghost, 1,         MPI_INT, send_proc, 606+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_S_inds, n_S_ghost, MPI_INT, send_proc, 707+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_S_ghost, Dim, rec_S_x,     send_proc, 808+send_to, MPI_COMM_WORLD ) ;

      // Add the particles to the current positions //
      insert_in_array( n_S_ghost, rec_S_inds, rec_S_x, x ) ;

      // Add particles to my_inds
      for ( i=ns_loc ; i < ns_loc + n_S_ghost ; i++ ) {
        my_inds[i] = rec_S_inds[i-ns_loc] ;
        local_flag[ rec_S_inds[i-ns_loc] ] = 1 ;
      }

      if ( mu != 0.0 ) {
        mpi_Recv_2d( n_S_ghost, Dim, rec_S_x,     send_proc, 828+send_to, MPI_COMM_WORLD ) ;
        insert_in_array( n_S_ghost, rec_S_inds, rec_S_x, mono_u ) ;
      }

      ns_loc += n_S_ghost ;

    }
  }//for ( send_proc=0:nprocs


  // Synchronize before moving to next loop
  MPI_Barrier( MPI_COMM_WORLD ) ;

  // Send the particles north //
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc + 1 ;
    if ( send_to >= nprocs ) send_to = 0 ;

    if ( send_proc == myrank ) {
      extract_from_array( Nbound_ct, send_N_inds, send_N_x, x ) ;
      for ( i=0 ; i<Nbound_ct ; i++ ) local_flag[ send_N_inds[i] ] = 0 ;

      MPI_Send( &Nbound_ct, 1, MPI_INT, send_to, 333+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_N_inds, Nbound_ct, MPI_INT, send_to, 444+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Nbound_ct, Dim, send_N_x, send_to, 555+send_to, MPI_COMM_WORLD ) ;
      
      // Send orientation vectors if necessary
      if ( mu != 0.0 ) {
        extract_from_array( Nbound_ct, send_N_inds, send_N_x, mono_u ) ;
        mpi_Send_2d( Nbound_ct, Dim, send_N_x, send_to, 528+send_to, MPI_COMM_WORLD ) ;
      }
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_N_ghost, 1,         MPI_INT, send_proc, 333+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_N_inds, n_N_ghost, MPI_INT, send_proc, 444+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_N_ghost, Dim, rec_N_x,     send_proc, 555+send_to, MPI_COMM_WORLD ) ;

      // Update positions positions //
      insert_in_array( n_N_ghost, rec_N_inds, rec_N_x, x ) ;

      for ( i=ns_loc ; i<ns_loc + n_N_ghost ; i++ ) {
        my_inds[i] = rec_N_inds[i-ns_loc] ;
        local_flag[ rec_N_inds[i-ns_loc] ] = 1 ;
      }
      
      if ( mu != 0.0 ) {
        mpi_Recv_2d( n_N_ghost, Dim, rec_N_x,     send_proc, 528+send_to, MPI_COMM_WORLD ) ;
        insert_in_array( n_N_ghost, rec_N_inds, rec_N_x, mono_u ) ;
      }     

      ns_loc += n_N_ghost ;
    }

  }//for ( send_proc=0:nprocs

  swap_t_out = time(0) ;
  swap_tot_time += swap_t_out - swap_t_in ;
  swap_partics_tot_time += swap_t_out - swap_t_in ;
}




// This sums the forces involving the ghost particles across
// the relevant processors. Occurs in two steps: 
// Forces are sent/received in the southern direction, then the
// northern direction
void communicate_ghost_forces() {

  if ( nprocs == 1 )
    return ; 


  swap_t_in = time(0);

  int i, send_proc, send_to, rec_from ;

  prepare_ghost_force_lists() ;

  // Send the forces south //
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc - 1 ;
    if ( send_to < 0 ) send_to = nprocs - 1 ;

    if ( send_proc == myrank ) {
      extract_from_array( Sbound_ct, send_S_inds, send_S_x, f ) ;

      MPI_Send( &Sbound_ct, 1, MPI_INT, send_to, 666+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_S_inds, Sbound_ct, MPI_INT, send_to, 777+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Sbound_ct, Dim, send_S_x, send_to, 888+send_to, MPI_COMM_WORLD ) ;
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_S_ghost, 1,         MPI_INT, send_proc, 666+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_S_inds, n_S_ghost, MPI_INT, send_proc, 777+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_S_ghost, Dim, rec_S_x,     send_proc, 888+send_to, MPI_COMM_WORLD ) ;

      // Add the forces to the current values //
      add_to_array( n_S_ghost, rec_S_inds, rec_S_x, f ) ;
    }
  }//for ( send_proc=0:nprocs


  // Syncronize before moving to nextx loop
  MPI_Barrier( MPI_COMM_WORLD ) ;


  // Send the forces north //
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc + 1 ;
    if ( send_to >= nprocs ) send_to = 0 ;

    if ( send_proc == myrank ) {
      extract_from_array( Nbound_ct, send_N_inds, send_N_x, f ) ;

      MPI_Send( &Nbound_ct, 1, MPI_INT, send_to, 101+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_N_inds, Nbound_ct, MPI_INT, send_to, 202+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Nbound_ct, Dim, send_N_x, send_to, 303+send_to, MPI_COMM_WORLD ) ;
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_N_ghost, 1,         MPI_INT, send_proc, 101+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_N_inds, n_N_ghost, MPI_INT, send_proc, 202+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_N_ghost, Dim, rec_N_x,     send_proc, 303+send_to, MPI_COMM_WORLD ) ;

      // Add the forces to the current values //
      add_to_array( n_N_ghost, rec_N_inds, rec_N_x, f ) ;
    }
  }//for ( send_proc=0:nprocs

  swap_t_out = time(0) ;
  swap_tot_time += swap_t_out - swap_t_in ;
  force_comm_tot_time += swap_t_out - swap_t_in ;

}





// This routine exchanges the indices and positions 
// of ghost particles. 
void swap_ghosts() {
  swap_t_in = time(0);

  if ( nprocs == 1 )
    return ;
  

  int i, j, rec_ct=0, send_to, rec_from, id, send_proc, rec_proc ;

  prepare_ghost_index_lists() ;


  // Send Southbound passengers //
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc - 1 ;
    if ( send_to < 0 ) send_to = nprocs - 1 ;

    if ( send_proc == myrank ) {
      MPI_Send( &Sbound_ct, 1,         MPI_INT, send_to, 100+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_S_inds,     Sbound_ct, MPI_INT, send_to, 200+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Sbound_ct, Dim, send_S_x,     send_to, 300+send_to, MPI_COMM_WORLD ) ;
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_S_ghost,  1,         MPI_INT, send_proc, 100+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_S_inds, n_S_ghost,  MPI_INT, send_proc, 200+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_S_ghost, Dim, rec_S_x,  send_proc, 300+send_to, MPI_COMM_WORLD ) ;
    }

  }
  
  // Set up first round of received particles //
  total_ghost = 0 ;
  for ( i=0 ; i<n_S_ghost ; i++ ) {
    id = rec_S_inds[i] ;
 
    if ( local_flag[id] == 1 ) continue ;

    for ( j=0 ; j<Dim ; j++ ) x[id][j] = rec_S_x[i][j] ;
    
    ghost_inds[ total_ghost ] = id ;
    total_ghost++ ;
  }




  // Set a barrier before sending/receiving N-bound passengers
  MPI_Barrier( MPI_COMM_WORLD ) ;


  // Send Northbound passengers
  for ( send_proc=0 ; send_proc < nprocs ; send_proc++ ) {

    send_to = send_proc + 1 ;
    if ( send_to >= nprocs ) send_to = 0 ;

    if ( send_proc == myrank ) {
      MPI_Send( &Nbound_ct, 1,         MPI_INT, send_to, 400+send_to, MPI_COMM_WORLD ) ;
      MPI_Send( send_N_inds,     Nbound_ct, MPI_INT, send_to, 500+send_to, MPI_COMM_WORLD ) ;
      mpi_Send_2d( Nbound_ct, Dim, send_N_x,     send_to, 600+send_to, MPI_COMM_WORLD ) ;
    }

    if ( send_to == myrank ) {
      MPI_Recv( &n_N_ghost,  1,        MPI_INT, send_proc, 400+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_N_inds, n_N_ghost, MPI_INT, send_proc, 500+send_to, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_N_ghost, Dim, rec_N_x, send_proc, 600+send_to, MPI_COMM_WORLD ) ;
    }

  }



  // Update final set of coordinates //
  for ( i=0 ; i<n_N_ghost ; i++ ) {
    id = rec_N_inds[i] ;
    
    if ( local_flag[id] == 1 ) continue ;

    for ( j=0 ; j<Dim ; j++ ) 
      x[id][j] = rec_N_x[i][j] ;

    ghost_inds[ total_ghost ] = id ;
    total_ghost++;
  }
  
  swap_t_out = time(0) ;
  swap_tot_time += swap_t_out - swap_t_in ;
  swap_ghosts_tot_time += swap_t_out - swap_t_in ;
}






// This routine updates all coordinates on the root process. It is
// useful for writing to the trajectory file, which is hopefully done
// at infrequent intervals.
void send_all_x_to_root() {

  if ( nprocs == 1 )
    return ;

  int i, send_proc ;
  for ( send_proc=1 ; send_proc < nprocs ; send_proc++ ) {

    if ( myrank == send_proc ) {
      extract_from_array( ns_loc, my_inds, send_S_x, x ) ;
      
      MPI_Send( &ns_loc, 1,      MPI_INT, 0, 1150+send_proc, MPI_COMM_WORLD ) ;
      MPI_Send( my_inds, ns_loc, MPI_INT, 0, 1250+send_proc, MPI_COMM_WORLD ) ;
      mpi_Send_2d( ns_loc, Dim, send_S_x, 0, 1350+send_proc, MPI_COMM_WORLD ) ;
    }


    else if ( myrank == 0 ) {
      MPI_Recv( &n_S_ghost, 1,         MPI_INT, send_proc, 1150+send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      MPI_Recv( rec_S_inds, n_S_ghost, MPI_INT, send_proc, 1250+send_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
      mpi_Recv_2d( n_S_ghost, Dim, rec_S_x, send_proc, 1350+send_proc, MPI_COMM_WORLD ) ;

      insert_in_array( n_S_ghost, rec_S_inds, rec_S_x, x ) ;
    
    }// myrank == 0

  }// for ( send_proc=1...

}

// Loops over the current procs ghost particles
// and builds the list of sites with non-zero forces
// to send to nearby procs. This routine currently
// sends all non-zero forces to all neighboring procs,
// which obv could be more efficient. PBCs annoying to
// deal with, and I'm just trying to debug at this point.
void prepare_ghost_force_lists() {
  int i,j,id ;
  double fsum  ;
  Sbound_ct = Nbound_ct = 0 ;

  for ( i=0 ; i<total_ghost ; i++ ) {
    id = ghost_inds[i] ;
    fsum = 0.0 ;
    for ( j=0 ; j<Dim ; j++ ) fsum += fabs(f[id][j]) ;

    if ( !local_flag[id] && fsum > 1.0E-6 ) {
      send_S_inds[ Sbound_ct ] = id ;
      Sbound_ct++ ;

      send_N_inds[ Nbound_ct ] = id ;
      Nbound_ct++ ;
    }
  }
}



// Find the indices for current processor to send both north and south //
// For a particle to be sent as a ghost particle it must be either:
// 1: within send_buff of the edge of the domain belonging to current proc
// or
// 2: Be bonded to a particle whose partner is not on the current proc
void prepare_ghost_index_lists() {

  int i, j, id, id2, k ;
  

  Sbound_ct = Nbound_ct = 0 ;

  if ( nprocs == 1 ) return ;
  
 
  int *N_flags, *S_flags ;
  N_flags = new int[ns_loc] ;
  S_flags = new int[ns_loc] ;


  for ( i=0 ; i<ns_loc ; i++ ) N_flags[i] = S_flags[i] = 0 ;


  // Tag particles within send_buff of boundary
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    
    if ( x[id][Dim-1] < z_min + send_buff ) 
      S_flags[i] = 1;

    if ( x[id][Dim-1] > z_max - send_buff ) 
      N_flags[i] = 1;
  }


  // Bonded particles that must be sent
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;

    for ( j=0 ; j<n_bonds[id]; j++ ) {
      id2 = bonded_to[id][j] ;

      // If the bond partner is on another proc, send "id"
      if ( local_flag[id2] == 0 ) continue ;

      S_flags[i] = 1 ;
      N_flags[i] = 1 ;
    }
  }

  for ( i=0 ; i<ns_loc ; i++ ) {
    int i1 = my_inds[i] ;

    for ( j=0 ; j<n_angles[i1] ; j++ ) {
      int jid[3], dd ;
      jid[0] = angle_first[i1][j] ;
      jid[1] = angle_mid[i1][j] ;
      jid[2] = angle_end[i1][j] ;

      if ( local_flag[jid[0]] && local_flag[jid[1]] && local_flag[jid[2]] ){ 
        continue ;
      }
      S_flags[i] = 1 ;
      N_flags[i] = 1 ;

    }// j=0:n_angles[id]

  }// i=0:ns_loc




  /////////////////////////////////////////
  // Build lists of all tagged particles //
  /////////////////////////////////////////
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;

    if ( S_flags[i] == 1 ) {
      send_S_inds[ Sbound_ct ] = id ;
      
      for ( j=0 ; j<Dim ; j++ ) send_S_x[ Sbound_ct ][j] = x[id][j] ;

      Sbound_ct += 1 ;
    }

    if ( N_flags[i] == 1 ) {
      send_N_inds[ Nbound_ct ] = id ;
      
      for ( j=0 ; j<Dim ; j++ ) send_N_x[ Nbound_ct ][j] = x[id][j] ;

      Nbound_ct += 1 ;
    }
  }// for i=0:ns_loc
  


  delete N_flags ;
  delete S_flags ;
}





// Takes values from f_short and adds them to values in f_long 
// at the memory locations given in inds
void add_to_array( int n, int *inds, double **f_short, double **f_long ) {
  int i, j, id ;
  for ( i=0 ; i<n ; i++ ) {
    id = inds[i] ;

    if ( local_flag[id] == 0 ) continue ;
    
    for ( j=0 ; j<Dim ; j++ ) f_long[id][j] += f_short[i][j] ;
  }

}


// Takes positions from x_short and places them in x_long at the memory 
// locations given in inds
void insert_in_array( int n, int *inds, double **x_short, double **x_long ) {
  int i, j, id ;
  for ( i=0 ; i<n ; i++ ) {
    id = inds[i] ;
    
    for ( j=0 ; j<Dim ; j++ ) x_long[id][j] = x_short[i][j] ;
  }

}



// Extracts the positions in inds from x_long and places them sequentially in x_short
void extract_from_array( int n, int *inds, double **x_short, double **x_long ) {
  int i, j, id ;
  for ( i=0 ; i<n ; i++ ) {
    id = inds[i] ;

    for ( j=0 ; j<Dim ; j++ ) x_short[i][j] = x_long[id][j] ;
  }
}


int is_partic_in_list( int n_list, int id, int *list ) {
  int i ;
  for ( i=0 ; i<n_list ; i++ ) {
    if ( list[i] == id ) return 1 ;
    }

  return 0;
}


// Finds the particles in the current processor's domain //
// This routine is only expected to be used upon         //
// initialization of the main code.
void find_my_particles() {

  int i ;

  ns_loc = 0 ;

  for ( i=0 ; i<nstot ; i++ ) {
    if ( x[i][Dim-1] > z_min && x[i][Dim-1] < z_max ) {
      my_inds[ ns_loc ] = i;
      local_flag[i] = 1 ;
      ns_loc++ ;
    }
    else
      local_flag[i] = 0 ;
  }
}




int remove_particles( int n_cur, int *cur_list, int n_kill, int *kill_list ) {

  int i, j, rm_id, k, n_removed = 0 ;

  for ( i=0 ; i<n_kill ; i++ ) {
    rm_id = kill_list[i] ;

    for ( j=0 ; j<n_cur ; j++ ) {
      if ( cur_list[j] == rm_id ) {

        for ( k=j ; k<n_cur ; k++ ) cur_list[k] = cur_list[k+1] ;
        
        n_removed++ ;
        break ;
      }
    }//for ( j=0 ; j<n_cur 

  }//for ( i=0:n_kill

  return n_cur - n_removed ;
}



