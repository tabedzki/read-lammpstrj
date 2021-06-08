#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
//#include "omp.h"
#include "mpi.h"
#include "fftw3-mpi.h" 

#define MPI_FLAG
#include "mpi_vars.h"

using namespace std ;




int mpi_Recv_2d( int ns, int dim, double **x, int origin, int tag, MPI_Comm comm ) {

  int i,j ;

  MPI_Recv( rec_xu, ns*dim, MPI_DOUBLE, origin, tag, comm, MPI_STATUS_IGNORE ) ;

  for ( i=0 ; i<ns ; i++ ) 
    for ( j=0 ; j<dim ; j++ ) 
      x[i][j] = rec_xu[i*dim + j] ;

  return 0 ;

}



int mpi_Send_2d( int ns, int dim, double **x, int dest, int tag, MPI_Comm comm ) {

  int i,j, ind=0 ;
  for ( i=0 ; i<ns ; i++ )
    for ( j=0 ; j<dim ; j++ )
      send_xu[ ind++ ] = x[i][j] ;

  MPI_Send( send_xu, ns*dim, MPI_DOUBLE, dest, tag, comm ) ;

  return 0 ;

}


// mpi_Bcast_posits: this routine is necessary because MPI expects to use contiguous
// memory when it sends/receives, which is screwed up with 2D arrays in C/C++. This 
// packs the data to be send into an array that is properly allocated, then sends it.
//
// Rob Riggleman   28 January 2018
//
// ns:   Number of sites to be send in current broadcast
// dim:  dimensionality of the system
// x:    ns*dim position array
// root: proc sending the data 
// comm: MPI_COMM_WORLD

int mpi_Bcast_posits( int ns, int dim, double **x, int root, MPI_Comm comm ) {

  int i,j, ind=0 ;
  
  for ( i=0 ; i<ns ; i++ )
    for ( j=0 ; j<dim ; j++ )
      tmp_x[ ind++ ] = x[i][j] ;

  MPI_Bcast( tmp_x, ns*dim, MPI_DOUBLE, root, comm ) ;

  for ( i=0 ; i<ns ; i++ ) 
    for ( j=0 ; j<dim ; j++ ) 
      x[i][j] = tmp_x[i*dim + j] ;

  return 0 ;

}




// d_ns: Dim * ns. Maximum size of a position array to send 
// Mtot: total grid points for possibly later needing to broadcast fields
int allocate_mpi_tmp( int d_ns, int Mtot ) {

  max_dof = d_ns ;
  tmp_x   = ( double* ) malloc( d_ns * sizeof(double) ) ;
  send_xu = ( double* ) malloc( d_ns * sizeof(double) ) ;
  rec_xu  = ( double* ) malloc( d_ns * sizeof(double) ) ;


  return ( 3*d_ns*sizeof(double) ) ;
}


