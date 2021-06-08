#include "globals.h"

void anneal_update( ) {
  int i ;

  cur_anneal_step++ ;

  ///////////////////////////////////////////////
  // Write out maximum of the structure factor //
  ///////////////////////////////////////////////
  double Imax = -23494.0,  kImax, k2, kv[Dim] ;
  char nm[40] ;

  for ( i=1 ; i<M ; i++ ) {
    if ( abs(avg_sk[0][i])/num_averages > Imax ) {
      Imax = abs(avg_sk[0][i]) ;
      kImax = sqrt( get_k( i , kv ) ) ;
    }
  }


  FILE *otp ;
  if ( cur_anneal_step == 1 ) 
    otp = fopen("chiAB_Imax_k.dat", "w" ) ;
  else
    otp = fopen("chiAB_Imax_k.dat", "a" ) ;

  fprintf(otp, "%lf %e %e\n", chiAB, Imax, kImax ) ;

  fclose( otp ) ;
  

  sprintf(nm, "avg_sk_cn%1.3lf.dat", chiAB ) ;
  
  for ( i=0 ; i<M ; i++ )
    ktmp2[i] = avg_sk[0][i] / num_averages ;
 
  write_kspace_data( nm , ktmp2 ) ;
  

  ///////////////////////////////
  // Zero the structure factor //
  ///////////////////////////////
  for ( i=0 ; i<M ; i++ ) 
    avg_sk[0][i] = 0.0 ;
  num_averages = 0.0 ;
  
  

  /////////////////
  // Update \chi //
  /////////////////
  if ( cur_anneal_step != n_anneal_pts ) {

    chiAB = anneal_chi[ cur_anneal_step ] ;

    fprintf(stderr, "chiAB changed to %lf!\n", chiAB) ;
    next_anneal_update = step + anneal_steps[ cur_anneal_step ] ;

  }


  if ( cur_anneal_step == 1 ) {
    Diff[1] = 0.1 ;
    //Diff[0] = Diff[1] = 0.1 ;
    printf("mobility update 1!\n");
  }

  else if ( cur_anneal_step == 2 ) {
    Diff[0] = 1.0 ;
    Diff[1] = 0.35 ;
    printf("mobility update 2!\n") ;
  }

}








void read_anneal( void ) {

  FILE *inp ;
  
  inp = fopen( "anneal.input", "r") ;
  if ( inp == NULL ) {
    do_anneal = 0 ;
    return ;
  }

  printf("Reading from anneal.input!!\nOnly affects \\chi_{AB}!\n\n") ;

  if ( !fscanf( inp , "%d\n", &n_anneal_pts ) )
    die("Failed to read first line of anneal.input");

  anneal_chi = ( double* ) calloc( n_anneal_pts , sizeof( double ) ) ;
  anneal_steps = ( int* ) calloc( n_anneal_pts, sizeof( int ) ) ;

  int i ;
  for ( i=0; i<n_anneal_pts ; i++ ) 
    if ( !fscanf( inp , "%lf %d\n", &anneal_chi[i] , &anneal_steps[i] ) )
      die("Failed to read line!\n") ;

  fclose(inp) ;

  do_anneal = 1 ;

  next_anneal_update = anneal_steps[0] ;
  chiAB = anneal_chi[0] ;

  cur_anneal_step = 0 ;


  nsteps = 0 ;
  for ( i=0 ; i<n_anneal_pts ; i++ ) 
    nsteps += anneal_steps[i] ;

  nsteps += 1 ;

  printf("\nmax steps updated to nsteps = %d!\n\n", nsteps ) ;

  sample_wait = 0 ;
}

