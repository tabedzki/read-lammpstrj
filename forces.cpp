#include <iomanip>
#include "globals.h"
#include "stiff-polymer-backbone.h"
void charge_grid( void ) ;
void bonds_gauss( void ) ;
void bonds_harmonic( void ) ;
void charge_forces( void ) ;
void communicate_ghost_forces( void ) ;
void angles( void ) ;
void angles_alt( void ) ;
void particle_orientations( void ) ;
void calc_S_conv_uG(void);
void ms_forces( void ) ;
void calc_scalar_parameters();



void forces() {

  int i,j, m, gind, t1, t2, id ;


  charge_grid() ;


  ///////////////////////////////////////////////
  // Reset the particle forces and grid grad w //
  ///////////////////////////////////////////////
  //Set forces for local particles to 0
  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
  //Set forces for ghost particles associated with current proc to 0
  for ( i=0 ; i<total_ghost ; i++ ) {
    id = ghost_inds[i] ;
    for ( j=0 ; j<Dim ; j++ )
      f[id][j] = 0.0 ;
  }
 
  //Sets gradients to 0 
  for ( i=0 ; i<ML ; i++ )
    for ( j=0 ; j<Dim ; j++ ) 
      gradwA[j][i] = gradwB[j][i] = gradwC[j][i] = gradwC[j][i] = gradwLC[j][i] = 0.0 ;



  //////////////////////////////////////////////////
  // Accumulate the monomer-monomer contributions //
  //////////////////////////////////////////////////
  
  // A acting on B, C //
  if (chiAC != 0.0 || chiAB != 0.0) {
    fftw_fwd(rho[0], ktmp);

    for (j = 0; j < Dim; j++) {
      for (i = 0; i < ML; i++) ktmp2[i] = grad_uG_hat[j][i] * ktmp[i];

      fftw_back(ktmp2, tmp);

      if (chiAB != 0.0)
        for (i = 0; i < ML; i++) {
          gradwB[j][i] += tmp[i] * chiAB / rho0;
        }
      if (chiAC != 0.0)
        for (i = 0; i < ML; i++) {
          gradwC[j][i] += tmp[i] * chiAC / rho0;
        }
    }
  }

  // B acting on A, C //
  if (chiAB != 0.0 || chiBC != 0.0) {
  fftw_fwd( rho[1] , ktmp ) ;
  
  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

      if ( chiAB != 0.0 )
    for ( i=0 ; i<ML ; i++ ) {
        gradwA[j][i] += tmp[i] * chiAB / rho0 ;
    }
      if ( chiBC != 0.0 )
    for ( i=0 ; i<ML ; i++ ) {
        gradwC[j][i] += tmp[i] * chiBC / rho0 ;
    }
  }}

  // C acting on A, B //
  if (chiAC != 0.0 || chiBC != 0.0) {
    fftw_fwd(rho[3], ktmp);

    for (j = 0; j < Dim; j++) {
      for (i = 0; i < ML; i++) ktmp2[i] = grad_uG_hat[j][i] * ktmp[i];

      fftw_back(ktmp2, tmp);

      if (chiAC != 0.0)
        for (i = 0; i < ML; i++) {
          gradwA[j][i] += tmp[i] * chiAC / rho0;
        }
      if (chiBC != 0.0)
        for (i = 0; i < ML; i++) {
          gradwB[j][i] += tmp[i] * chiBC / rho0;
        }
    }
  }

  // Compressibility contribution //
  if (do_field == 1){
  for ( i=0 ; i<ML ; i++ )
    tmp[i] = rhot[i] ;

  fftw_fwd( tmp , ktmp ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    for ( i=0 ; i<ML ; i++ )
      ktmp2[i] = grad_uG_hat[j][i] * ktmp[i] ;

    fftw_back( ktmp2 , tmp ) ;

    for ( i=0 ; i<ML ; i++ ) {
      gradwA[j][i]  += tmp[i] * kappa / rho0 ;
      gradwB[j][i]  += tmp[i] * kappa / rho0 ;
      gradwC[j][i]  += tmp[i] * kappa / rho0 ;
      gradwLC[j][i] += tmp[i] * kappa / rho0 ;
    }

  }
  }

  ///////////////////////////////////////////
  // Accumulate the particle contributions //
  ///////////////////////////////////////////
  if ( nP > 0 ) {

    // Particle-particle //
    fftw_fwd( rho[2] , ktmp ) ;

    for ( j=0 ; j<Dim ; j++ ) {
      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uP_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;
    }

    // Particles acting on monomers //
    for ( j=0 ; j<Dim ; j++ ) {
      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ ) {
        gradwA[j][i]  += tmp[i] * kappa / rho0 ;
        gradwLC[j][i] += tmp[i] * kappa / rho0 ;
        gradwC[j][i]  += tmp[i] * ( kappa + ( A_partics ? chiAC : 0.0 ) ) / rho0 ;
        gradwB[j][i]  += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;
      }
    }

    // A MLonomers acting on particles //
    fftw_fwd( rho[0] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++ ) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * kappa / rho0 ;

    }

    // B MLonomers acting on particles //
    fftw_fwd( rho[1] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAB : 0.0 ) ) / rho0 ;

    }


    // C MLonomers acting on particles //
    fftw_fwd( rho[3] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa + ( A_partics ? chiAC : 0.0 ) ) / rho0 ;

    }


    // LC chains acting on particles //
    fftw_fwd( rho[4] , ktmp ) ;
    for ( j=0 ; j<Dim ; j++) {

      for ( i=0 ; i<ML ; i++ )
        ktmp2[i] = grad_uPG_hat[j][i] * ktmp[i] ;

      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<ML ; i++ )
        gradwP[j][i] += tmp[i] * ( kappa ) / rho0 ;

    }
  }// if ( nP > 0 )



      

  //////////////////////////////////////////////////////
  // Accumulate the nonbonded forces on each particle //
  //////////////////////////////////////////////////////
  for ( i=0 ; i<ns_loc+total_ghost ; i++ ) {

    if ( i < ns_loc )
      id = my_inds[i] ;
    else
      id = ghost_inds[ i-ns_loc ] ;


    for ( m=0 ; m < grid_per_partic ; m++ ) {

      // Note that in charge_grid(), grid_inds is populated with 
      // indices in the range [0,ML]. 
      // if gind == -1, then grid location is skipped.
      gind = grid_inds[ id ][ m ] ;


      if ( gind == -1 )
        continue ;


      for ( j=0 ; j<Dim ; j++ ) {
        if ( mol_type[id] == 0 )
          f[id][j] -= gradwA[ j ][ gind ] * grid_W[id][m] ;
        else if ( mol_type[id] == 1 )
          f[id][j] -= gradwB[ j ][ gind ] * grid_W[id][m] ;
        else if ( mol_type[id] == 2 )
          f[id][j] -= gradwP[ j ][ gind ] * grid_W[id][m] ;
        else if ( mol_type[id] == 3 )
          f[id][j] -= gradwC[ j ][ gind ] * grid_W[id][m] ;
        else if ( mol_type[id] == 4 )
          f[id][j] -= gradwLC[ j ][ gind ] * grid_W[id][m] ;
      }
    }

    for ( j=0 ; j<Dim ; j++ )
      f[id][j] *= grid_vol ;
  }


  ////////////////////////////
  // Call the bonded forces //
  ////////////////////////////
  bonds_gauss();
  // angles();
  angles_alt();

  //////////////////////////////////
  // Shear-able stretchable model //
  //////////////////////////////////

  if (do_dsswlc == 1){
  stiff::calculate();
  }

  if ( mu != 0.0 || print_vectors ==1 ) {
    calc_S_conv_uG() ;
    ms_forces() ;
  }


  if ( nprocs > 1 )
    communicate_ghost_forces() ;

}
