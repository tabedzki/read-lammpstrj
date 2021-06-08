#include "globals.h"
#include "dump.h"
void lagrange_get_weights( double , double , double* ) ;
void spline_get_weights( double , double , double* ) ;
void add_segment( int ) ;
  double Lloc[3];

extern int t;
double pbc_mdr2( double *r1 , double *r2 , double *dr ) {

  Lloc[0] = L[0][0];
  Lloc[1] = L[0][1];
  Lloc[2] = L[0][2];
  int j ;
  double mdr2 = 0.0 ;

  for ( j=0 ; j<Dim ; j++ ) {
    dr[j] = r1[j] - r2[j] ;
    
    if ( dr[j] > Lh[j] )
      dr[j] -= Lloc[j] ;

    else if ( dr[j] <= -Lh[j] )
      dr[j] += Lloc[j] ;

    mdr2 += dr[j] * dr[j] ;
  }

  return mdr2 ;

}

// Random vector on the unit sphere/circle //
void random_u( double u[Dim] ) {

  int i ;
  double mdr2 = 0.0 , mdr ;
  while ( mdr2 == 0.0 ) {

    for ( i=0 ; i<Dim ; i++ ) {
      u[i] = 2.0 * ( ran2() - 0.5 ) ;
      mdr2 += u[i] * u[i] ;
    }

    if ( mdr2 > 1.0 )
      mdr2 = 0.0 ;
  }

  mdr = sqrt( mdr2 ) ;

  for ( i=0 ; i<Dim ; i++ )
    u[i] = u[i] / mdr ;

}

void convolve_fields( double *in1 , double *in2 , double *out ) {
  int i ;

  fftw_fwd( in1 , ktmp ) ;
  fftw_fwd( in2 , ktmp2 ) ;

  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , out ) ;

}

void field_gradient( double *inp , double *out , int dir ) {

  int i ;
  double k2, kv[Dim] ;
  fftw_fwd( inp , ktmp ) ;

  for ( i=0 ; i<M ; i++ ) {
    k2 = get_k_alias( i , kv ) ;

    ktmp[i] *= I * kv[ dir ] ;
  }

  fftw_back( ktmp , out ) ;

}

void field_gradient_cdif( double *inp , double *out , int dir ) {

  
  int i, j, nx2, nx, px, px2, nn[Dim], ntmp[Dim] ;
  
  for ( i=0 ; i<M ; i++ ) {

    unstack( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      ntmp[j] = nn[j] ;

    
    ntmp[ dir ] = nn[ dir ] - 2 ;
    if ( ntmp[ dir ] < 0 ) ntmp[ dir ] += Nx[ dir ] ;
    nx2 = stack( ntmp ) ;

    ntmp[ dir ] = nn[ dir ] - 1 ;
    if ( ntmp[ dir ] < 0 ) ntmp[ dir ] += Nx[ dir ] ;
    nx = stack( ntmp ) ;


    ntmp[ dir ] = nn[ dir ] + 1 ;
    if ( ntmp[ dir ] >= Nx[ dir ] ) ntmp[ dir ] -= Nx[ dir ] ;
    px = stack( ntmp ) ;

    ntmp[ dir ] = nn[ dir ] + 2 ;
    if ( ntmp[ dir ] >= Nx[ dir ] ) ntmp[ dir ] -= Nx[ dir ] ;
    px2 = stack( ntmp ) ;

    out[i] =  ( inp[ nx2 ] - 8.0 * inp[ nx ] + 8.0 * inp[ px ] - inp[ px2 ] ) 
      / ( 12.0 * dx[ dir ] ) ;

  }
 

}



int stack_input(int x[Dim], int Nxx[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nxx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nxx[1])*Nxx[0] );
}

// Stacks x using only local values
int stack_local(int x[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}

// Stacks vector x into 1D array index in [ 0, M ]
int stack( int x[Dim] ) { 
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}

void unstack_input(int id, int nn[Dim], int Nxx[Dim]) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nxx[0];
    nn[0] = (id - nn[1]*Nxx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nxx[1]/Nxx[0];
    nn[1] = id/Nxx[0] - nn[2]*Nxx[1];
    nn[0] = id - (nn[1] + nn[2]*Nxx[1])*Nxx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

// Receives index id in [0 , M ] and makes array
// nn[Dim] in [ 0 , Nx[Dim] ]
void unstack(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nx[0];
    nn[0] = (id - nn[1]*Nx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nx[1]/Nx[0];
    nn[1] = id/Nx[0] - nn[2]*Nx[1];
    nn[0] = id - (nn[1] + nn[2]*Nx[1])*Nx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

void get_r( int id , double r[Dim] ) {
  int i, id2, n[Dim];

  unstack(id, n);


  for ( i=0; i<Dim; i++) {
    r[i] = dx[i] * double( n[i] );
    
    if ( r[i] > Lloc[i]/2.0 )
      r[i] -= Lloc[i];
    else if ( r[i] <= -Lloc[i]/2.0 )
      r[i] += Lloc[i];
 
  }
}

double get_k_alias( int id , double k[Dim] ) {

  double kmag = 0.0;
  int i, id2, n[Dim] , j , has_nyquist = 0;
  for ( i=0 ; i<Dim ; i++ )
    if ( Nx[i] % 2 == 0 )
      has_nyquist = 1;

  unstack(id, n);

  if ( Nx[0] % 2 == 0 && n[0] == Nx[0] / 2 )
    k[0] = 0.0 ;
  else if ( double(n[0]) < double(Nx[0]) / 2.)
   k[0] = 2*PI*double(n[0])/Lloc[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/Lloc[0];

  if (Dim>1) {
    if ( Nx[1] % 2 == 0 && n[1] == Nx[1] / 2 )
      k[1] = 0.0 ;
    else if ( double(n[1]) < double(Nx[1]) / 2.)
      k[1] = 2*PI*double(n[1])/Lloc[1];
    else
      k[1] = 2*PI*double(n[1]-Nx[1])/Lloc[1];
  }

  if (Dim==3) {
    if ( Nx[2] % 2 == 0 && n[2] == Nx[2] / 2 )
      k[2] = 0.0 ;
    else if ( double(n[2]) < double(Nx[2]) / 2.)
      k[2] = 2*PI*double(n[2])/Lloc[2];
    else
      k[2] = 2*PI*double(n[2]-Nx[2])/Lloc[2];
  }

  // Kills off the Nyquist modes
  if ( id2 != 0 && has_nyquist ) {
    for ( i=0 ; i<Dim ; i++ ) {
      if ( k[i] == 0.0 ) {
        for ( j=0 ; j<Dim ; j++ )
          k[j] = 0.0 ;
        kmag = 0.0;
        break;
      }
    }
  }
  
  for (i=0; i<Dim; i++)
    kmag += k[i]*k[i];

  return kmag;

}



// Receives index id in [ 0 , ML ] and returns 
// proper k-value, whether running in parallel or not
double get_k(int id, double k[Dim]) {

  double kmag = 0.0;
  int i, id2, n[Dim];

  unstack(id, n);

  for ( i=0 ; i<Dim ; i++ ) {
    if ( double( n[i] ) < double( Nx[i] ) / 2. )
      k[i] = 2 * PI * double( n[i] ) / L[t][i] ;
    else
      k[i] = 2 * PI * double( n[i] - Nx[i] ) / L[t][i] ;

    kmag += k[i] * k[i] ;
  }

  return kmag;

}


// Sets the average of tp to zero
void zero_average(double* tp) {

  int i;

  double integ;
  
  integ = integrate(tp);

  integ *= (1.0 / V);

  for (i=0; i<M; i++)
    tp[i] -= integ;

}

