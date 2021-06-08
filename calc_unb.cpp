#include "globals.h"
double calc_ms_energy(void) ;


void calc_Unb() {

  int i ;
  // Take Gaussian potential to k-space
  fftw_fwd( uG , ktmp2 ) ;
  
  
  // Polymer kappa contribution //
  for ( i=0 ; i<ML ; i++ ) tmp2[i] = rhot[i] - rho0 ;

  fftw_fwd( tmp2 , ktmp ) ;

  for ( i=0 ; i<ML ; i++ ) ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

  for ( i=0 ; i<ML ; i++ ) tmp[i] *= tmp2[i] ;

  Ukappa = integrate( tmp ) * kappa / 2.0 / rho0 ;


  if ( mu != 0.0 )
    U_ms = calc_ms_energy() ;
  else
    U_ms = 0.0 ;

}

double calc_ms_energy() {

  int i,j,m,n ;
  double ddot ;

  for (i = 0; i < ML; i++) {
      tmp[i] = 0.0;
      for (m = 0; m < Dim; m++) {
          for (n = 0; n < Dim; n++){
              tmp[i] += S_field[i][m][n] * S_conv_u[i][n][m];
          }
      }
  }

  return ( integrate(tmp) * -mu / 2.0 / rho0 ) ;

}

