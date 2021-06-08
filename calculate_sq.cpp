#include <stdio.h>
#include <complex>
#include "stdlib.h"
#include "math.h"
#include "dump.h"
#include "globals.h"

using namespace std ;
void field_gradient( double *inp , double *out , int dir ) {

  int i ;
  double k2, kv[Dim] ;
  fftw_fwd( inp , ktmp ) ;

  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k_alias( i , kv ) ;

    ktmp[i] *= I * kv[ dir ] ;
  }

  fftw_back( ktmp , out ) ;

}

