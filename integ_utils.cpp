#include "globals.h"

double integrate( double *inp ) {
  double sum = 0.0 ;
  int i ;

#pragma omp parallel for reduction(+:sum)
  for ( i=0 ; i<M ; i++ )
    sum += inp[i] ;

  for ( i=0 ; i<Dim ; i++ )
    sum *= dx[i] ;

  return sum ;
}
