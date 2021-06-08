#define PI 3.141592653589793238462643383
#define MAIN
#include <stdio.h>
#include <complex>
#include "stdlib.h"
#include "math.h"
#include "dump.h"
#include "globals.h"
using namespace std ;

void qs(double*,int*,int,int);
void fft_init();

int read_dump_traj(char*, int , int);
void charge_grid();
void log_space(int,int,int);
double legendre(double);
double integ(int, double, double*);
int t;

int main(int argc,char *argv[])
{

  if (argc < 7) {
    printf("Usage: ./a.out [input.lammpstrj] [frmin] [frmax] [Nx] [Ny] [Nz] \n");
    exit(1);
  }


  if ( argc == 7 && !strcmp( "-nt" , argv[1] ) ) {
    nthreads = atoi( argv[2] ) ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }
  else {
    nthreads = 1 ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }

  int MAXFRAMES = 2500, MAXBIN = 2500, nrbins;
  int i,m,j,k,frs,delt,pflag=0 ;
  double **bx, dr, drbin = 0.002, sf=1.0,dy,dz;
  char name[20];  
  
  I = complex<double>( 0.0 , 1.0 ) ;

  static int n=1;
  int frmin=0, frmax ;
  frmin = atoi(argv[2]);
  frmax = atoi(argv[3]);
  Nx[0] = atoi(argv[4]);
  Nx[1] = atoi(argv[5]);
  Nx[2] = atoi(argv[6]);

  int ML = Nx[0] * Nx[1] * Nx[2];
  M=ML;
 
 ntypes=5;
 std::cout<<"I allocate this much: "<< ML<<endl;

  frs = frmax - frmin;
  // grid_per_partic  = 3;
  pmeorder = 3;

   rhot = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhoha = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhohb = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhohc = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhoda = ( double* ) calloc( ML , sizeof( double ) ) ;
  rhodb = ( double* ) calloc( ML , sizeof( double ) ) ;
  complex<double> *ktmp  = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
  complex<double> *ktmp2  = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
  complex<double> *ktmp3  = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
  complex<double> *ktmp_avg_a  = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
  complex<double> *ktmp_avg_b  = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;

  rho = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  rho_hat = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  avg_sk = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  w = ( double** ) calloc( ntypes , sizeof( double* ) ) ;

  for ( i=0 ; i<ntypes ; i++ ) {
    rho_hat[i] = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
    rho[i] = ( double * ) calloc( ML , sizeof( double ) ) ;
    avg_sk[i] = ( complex<double> * ) calloc( ML , sizeof( complex<double> ) ) ;
    w[i] = ( double* ) calloc( ML , sizeof( double ) ) ;
  }

  rhogb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoga_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoha_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhohb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoda_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhodb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhop_t =  ( double** ) calloc( nthreads , sizeof( double* ) ) ;

  cout<<ML<<" "<<M<<endl;
  for ( i=0 ; i<nthreads ; i++ ) {
    rhogb_t[i]  = ( double* ) calloc(M , sizeof( double ) ) ;
    rhoga_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhoha_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhohb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhoda_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhodb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhop_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
  }


  
  cout<<rhoda[0] << endl;
  rhoda[0]=0;
  cout<<rhodb[0] << endl;

  frs = read_dump_traj(argv[1], frmin, frmax);

  grid_inds = ( int** ) calloc( ns , sizeof( int* ) );
  grid_W = ( double** ) calloc( ns , sizeof( double* ) ) ;

  for ( i=0 ; i<ns ; i++ ) {
    grid_inds[i] = ( int* ) calloc( grid_per_partic , sizeof( int ) ) ;
    grid_W[i] = ( double* ) calloc( grid_per_partic , sizeof( double ) ) ;
  }
  
  nstot=ns;
  printf("%d frames\n",frs);
  bx = (double**) calloc(frs,sizeof(double*));
  for (t=0; t<frs; t++) {
    bx[t] = (double*) calloc(3,sizeof(double));
    bx[t][0] = L[t][0] ;
    bx[t][1] = L[t][1] ;
    bx[t][2] = L[t][2] ;

  M=1;
  V=1;
  cout<<M<<endl;
  for ( j=0 ; j<Dim ; j++ ) {
    V *= L[t][j] ;
  // cout<<2<<endl;
    M *= Nx[j] ;
  // cout<<3<<endl;
    // ML *= Nx[t][j] ;
    dx[j] = L[t][j] / double( Nx[j] ) ;
  // cout<<4<<endl;
    Lh[j] = 0.5 * L[t][j] ;
  // cout<<5<<endl;
    grid_per_partic *= ( pmeorder + 1 ) ;
  // cout<<6<<endl;
      /* printf("dx: %lf\n" , dx[j] ) ; */
  }
  gvol = V / double( M ) ;
  }
  ML=M;

  // Reconnect the trajectories //
  // for (t=1; t<frs; t++) {
  //   for (i=0; i<ns ; i++) {
  //     for (j=0; j<3; j++) {
  //       double bl = 0.5*( bx[t][j] + bx[t-1][j] );

  //       dx = (x[t][i][j] - x[t-1][i][j]) ;
  //       while ( dx > bl/2. ) {
  //         x[t][i][j] -= bl ;
  //         dx = (x[t][i][j] - x[t-1][i][j]) ;
  //       }
  //       while ( dx <= -bl/2. ) {
  //         x[t][i][j] += bl ;
  //         dx = (x[t][i][j] - x[t-1][i][j]) ;
  //       }

  //       dx = (x[t][i][j] - x[t-1][i][j]) ;
  //       if ( dx > 0.5 * bl || dx < -0.5*bl ) {
  //         printf("You messed up! %lf %d %d %d\n" , dx , t,k,j);
  //         exit(1);
  //       }
 
  //     }
  //   }
  // }

  printf("Particles' trajectories connected in time\n"); fflush( stdout );
  rhodb[0]=0.0;
  for (t=0; t<frs; t++) {
    printf("\rTimestep %d out of %d", t, frs); fflush( stdout );
    // if (t==0)
    // cout<<"first charge"<<endl;
    // if (t==1)
    // cout<<"second charge"<<endl;
    charge_grid();

    // printf("yellow\n"); fflush( stdout );
    fft_init();
    // printf("green\n"); fflush( stdout );
    // std::cout<<"rho"<<*rho[0] <<std::endl;
    // std::cout<<"rho"<<rho[0][0] <<std::endl;
    // for (i=0;i<ML;i++){
    //   std::cout<<rho[1][i]<<std::endl;
    // }
    // std::cout<<"rho"<<ktmp[0] <<std::endl;
    // write_grid_data("rho.dat", rho[0]);
          fftw_fwd( rho[0], ktmp ) ;

    // printf("brown\n"); fflush( stdout );
        for ( i=0 ; i<M ; i++ ){
          ktmp3[i] = ktmp[i] * conj(ktmp[i]) ;
          ktmp_avg_a[i] += ktmp3[i];
        }
        
        //  write_kspace_data( "a_conv_a.dat", ktmp3 ) ;

          fftw_fwd( rho[1], ktmp ) ;

    // printf("brown\n"); fflush( stdout );
        for ( i=0 ; i<M ; i++ ){
          ktmp3[i] = ktmp[i] * conj(ktmp[i]) ;
          ktmp_avg_b[i] += ktmp3[i];
        }
        
        //  write_kspace_data( "b_conv_b.dat", ktmp3 ) ;
        //   fftw_fwd( rho[1], ktmp2 ) ;

        // for ( i=0 ; i<M ; i++ )
        //   ktmp3[i] = ktmp2[i] * conj(ktmp2[i]) ;
        
        //  write_kspace_data( "sk", ktmp3 ) ;


  }


  t--;

        for ( i=0 ; i<M ; i++ ){
          ktmp_avg_a[i] /= frs;
          ktmp_avg_b[i] /= frs;
        }

        write_kspace_data( "a_conv_a.avg", ktmp_avg_a);
        write_kspace_data( "b_conv_b.avg", ktmp_avg_b);

  return 0;
}

