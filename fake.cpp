#define XYZ
#include "dump.h"
#include "globals.h"
#include "stdlib.h"
#include "stdio.h"
#include <limits.h>
#include <float.h>

int read_dump_traj(char *nm, int frmin , int frmax) {

  int i, j, frs, di, ind ;
  char tt[80];

  FILE *inp;
  inp = fopen(nm,"r");
  if (inp==NULL) { 
    printf("Failed to open %s!\n", nm);
    exit(1);
  }
  int MAX = frmax - frmin, t0, t1  ;
  x = (double***) calloc(MAX,sizeof(double**));
  L = ( double** ) calloc( MAX , sizeof( double* ) ) ;
  frs = MAX;

  fgets( tt , 80 , inp ) ;
  fscanf(inp,"%d\n",&t0);
  fgets( tt , 80 , inp ) ;
  fscanf(inp,"%d\n",&ns);
  for ( i=0 ; i<6+ns ; i++ )
    fgets(tt,80,inp) ;
  fscanf( inp,"%d\n", &t1) ;
  save_freq = t1 - t0 ;
  printf("Figured saving frequency to be %d\n", save_freq ) ;
  fclose(inp);

  // double frmin_index = (double) DBL_MAX;
  // double frmax_index = (double) -DBL_MAX;
  // double tmp = (double) DBL_MAX;
  // inp = fopen(nm,"r");
  //   for ( j=0 ; j < 9 ; j++ )
  //     fgets( tt , 80 , inp );
  //   for ( j=0 ; j < ns ; j++ ){
  //   fscanf( inp, "%lf" , &tmp ) ;
  //   frmin_index = min(tmp, frmin_index);
  //   frmax_index = max(tmp, frmax_index);
  //     fgets( tt , 80 , inp );
  //     }
  //   fclose(inp);

  inp = fopen( nm, "r" );
  for ( i=0 ; i < frmin ; i++ ) {
    for ( j=0 ; j < ns+9 ; j++ )
      fgets( tt , 80 , inp );
  }

  for ( j=0 ; j<MAX ; j++ ) {
    fgets(tt,80,inp);
    fgets(tt,80,inp);
    fgets(tt,80,inp);
    fgets(tt,80,inp);
    fgets(tt,80,inp);

    // Allocate this frame's memory //

    L[j] = ( double* ) calloc( 3 , sizeof( double ) ) ;
    // x[j] = (double**) calloc(frmax_index,sizeof(double*));
    x[j] = (double**) calloc(ns,sizeof(double*));
    // cout<<frmax_index<<endl;
    // for (i=0; i<frmax_index; i++)
    for (i=0; i<ns; i++)
      x[j][i] = (double*) calloc(3,sizeof(double));


    // Allocate memory for the particle types //
    if (j==0) {
      tp = (int*) calloc(ns,sizeof(int));
      // for (i=0; i<ns; i++)
      //   tp[i] = (int) calloc(1,sizeof(int));
    }

    double xlo, ylo, zlo, xhi, yhi, zhi;
    fscanf( inp, "%lf %lf\n" , &xlo, &xhi ) ;
    L[j][0] = xhi - xlo ;
    fscanf( inp, "%lf %lf\n" , &ylo, &yhi ) ;
    L[j][1] = yhi - ylo ;
    fscanf( inp, "%lf %lf\n" , &zlo, &zhi ) ;
    L[j][2] = zhi - zlo ;

    // std::cout<<L[j][0]<<std::endl;
    // std::cout<<L[j][1]<<std::endl;
    // std::cout<<L[j][2]<<std::endl;

    fgets(tt,80,inp);

    // Read the frame
    for (i=0; i<ns; i++) {
      fscanf(inp,"%d",&ind ) ;
      ind = ind - 1;
      // ind = ind - frmin_index ;
      fscanf(inp,"%d",&tp[ind]) ;
      fscanf(inp,"%d",&di ) ;

      fscanf(inp,"%lf",&x[j][ind][0]);if (feof(inp)) break;
      fscanf(inp,"%lf",&x[j][ind][1]);if (feof(inp)) break;
      fscanf(inp,"%lf",&x[j][ind][2]);if (feof(inp)) break;
      fgets(tt,80,inp);if (feof(inp)) break;
      x[j][ind][0] -= xlo;
      x[j][ind][1] -= ylo;
      x[j][ind][2] -= zlo;
    }

    printf("\rframe %d of %d MAX read", j, MAX); fflush(stdout);
  }
  printf("\n");
  printf("Read %d frames with %d sites\n",frs,ns);
  fclose(inp);
  int id=9;
  j=0;
  int t=0;
  cout<<x[t][id][j]<<'\t'<<dx[j]<<endl;
    // die("bb");
  return frs;

}
