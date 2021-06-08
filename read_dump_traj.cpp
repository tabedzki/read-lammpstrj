#define XYZ
#include "dump.h"
#include "globals.h"
#include "stdlib.h"
#include "stdio.h"
#include <limits.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

unsigned int FileRead( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

unsigned int CountLines( const vector <char> & buff, int sz ) {
    int newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
        if ( p[i] == '\n' ) {
            newlines++;
        }
    }
    return newlines;
}


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

  int count=0;
  cout << "buffer\n";
        const int SZ = 1024 * 1024;
        std::vector <char> buff( SZ );
        ifstream ifs( nm );
        while( int cc = FileRead( ifs, buff ) ) {
            count += CountLines( buff, cc );
        }

  ifs.close();
  // char c;
  // inp = fopen(nm,"r");
  //   for (c = getc(inp); c != EOF; c = getc(inp)) 
  //     if (c == '\n') // Increment count if this character is newline 
  //         count = count + 1; 

  int total_timesteps = count / (ns+9);

  if (frmax > total_timesteps){
    printf("frmax exceeds total timesteps.\nOverwriting frmax to the frame %d.\n",
      total_timesteps);
      frmax = total_timesteps;
    MAX=total_timesteps-frmin;
    frs=total_timesteps-frmin;
  }
  cout<<"MAX: "<<MAX<<"\tfrs: "<<frs<<"\ttotal_timesteps "<<total_timesteps<<endl;
  cout<<"frmin: "<<frmin<<"\tfrmax: "<<frmax<<endl;

  if (frs<0){
    frs = -frs;
    printf("Specified number of frames is negative. Taking the last %d frames\n",
         frs);
    frmin = total_timesteps - frs;
    frmax = total_timesteps;
    MAX=frs;
  }
  cout<<"frmin: "<<frmin<<"\tfrmax: "<<frmax<<endl;
  

  double frmin_index = (double) DBL_MAX;
  double frmax_index = (double) -DBL_MAX;
  double tmp = (double) DBL_MAX;
  inp = fopen(nm,"r");
    for ( j=0 ; j < 9 ; j++ )
      fgets( tt , 80 , inp );
    for ( j=0 ; j < ns ; j++ ){
    fscanf( inp, "%lf" , &tmp ) ;
    frmin_index = min(tmp, frmin_index);
    frmax_index = max(tmp, frmax_index);
      fgets( tt , 80 , inp );
      }
    fclose(inp);

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
    x[j] = (double**) calloc(frmax_index,sizeof(double*));
    // x[j] = (double**) calloc(ns,sizeof(double*));
    // cout<<frmax_index<<endl;
    for (i=0; i<frmax_index; i++)
    // for (i=0; i<ns; i++)
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

    // double xdiff = xhi-xlo;
    // double ydiff = yhi-ylo;

    // double zdiff = zhi-zlo;

    // double diff[3];
    // diff[0]=xdiff;
    // diff[1]=ydiff;
    // diff[2]=zdiff;
    // cout << diff[0] << diff[1]<< diff[2]<<endl;
    
    // Read the frame
    for (i=0; i<ns; i++) {
      fscanf(inp,"%d",&ind ) ;
      // ind = ind - 1;
      ind = ind - frmin_index ;
      fscanf(inp,"%d",&tp[ind]) ;
      fscanf(inp,"%d",&di ) ;

      fscanf(inp,"%lf",&x[j][ind][0]);if (feof(inp)) break;
      fscanf(inp,"%lf",&x[j][ind][1]);if (feof(inp)) break;
      fscanf(inp,"%lf",&x[j][ind][2]);if (feof(inp)) break;
      fgets(tt,80,inp);if (feof(inp)) break;
      x[j][ind][0] -= xlo;
      x[j][ind][1] -= ylo;
      x[j][ind][2] -= zlo;
      // for ( int llj=0 ; llj<Dim ; llj++ ) { 
	    //   if ( x[j][ind][llj] > diff[llj]  ) 
		  //     x[j][ind][llj] -= diff[llj];

	    //   else if ( x[j][ind][j] < 0.0 )
		  //     x[j][ind][llj] += diff[llj] ;
      // }   

    }

    printf("\rframe %d of %d MAX read", j, MAX); fflush(stdout);
  }
  printf("\n");
  printf("Read %d frames with %d sites\n",frs,ns);
  // fclose(inp);
  return frs;

}
