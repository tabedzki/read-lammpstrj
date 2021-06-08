#include "globals.h"

#include <sys/stat.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
// void write_stress( ) {

//   int i, j, k ;

//   FILE *otp ;

//   if ( step == 0 )
//     otp = fopen( "stress.dat" , "w" ) ;
//   else
//     otp = fopen( "stress.dat" , "a" ) ;
// /*  if ( step <= print_freq )
//     otp1 = fopen( "ete.dat" , "w" ) ;
//   else
//     otp1 = fopen( "ete.dat" , "a" ) ;
// */

// //  for ( i=0 ; i<buff_ind ; i++ ) {
//     // Write diagonals first //



//     for ( j=0 ; j<Dim ; j++ )
//         fprintf( otp , "%lf " , Stress_nb[j][j] ) ;

//   fprintf( otp , "%lf " , Stress_nb[0][1] ) ;


//     for ( j=0 ; j<Dim ; j++ )
//       fprintf( otp , "%lf " , Stress_bonds[j][j] ) ;

//         fprintf( otp , "%lf " , Stress_bonds[0][1] ) ;
//     //fprintf( otp , "%lf " , cum_pr*double(step - pre_equil_steps)/double(stress_freq));
    
//  //   fprintf( otp1 , "%lf " , sts_buf[i][Dim][Dim] ) ;
//     fprintf( otp , "\n" ) ;
//  //   fprintf( otp1 , "\n" ) ;
// //  }

//   fclose( otp ) ;
//  // fclose( otp1 ) ;
//  // buff_ind = 0 ;

// }
/*void write_stress( ) {

  int i, j, k ;

  FILE *otp,*otp_pp,*otp_ng ;

  if ( step <= print_freq )
    otp = fopen( "stress_euler.dat" , "w" ) ;
  else
    otp = fopen( "stress_euler.dat" , "a" ) ;

  if ( step <= print_freq )
    otp_pp = fopen( "stress_PP.dat" , "w" ) ;
  else
    otp_pp = fopen( "stress_PP.dat" , "a" ) ;

  if ( step <= print_freq )
    otp_ng = fopen( "stress_ga.dat" , "w" ) ;
  else
    otp_ng = fopen( "stress_ga.dat" , "a" ) ;



  for ( i=0 ; i<buff_ind ; i++ ) {
    // Write diagonals first //
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , sts_buf[i][j][j] ) ;

     for ( j=0 ; j<Dim ; j++ )
      fprintf( otp_pp , "%lf " , sts_buf_pp[i][j][j] ) ;

      for ( j=0 ; j<Dim ; j++ )
      fprintf( otp_ng , "%lf " , sts_buf_ng[i][j][j] ) ;

 

    for ( j=0 ; j<Dim  ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
        fprintf( otp , "%lf " , sts_buf[i][j][k] ) ;

    for ( j=0 ; j<Dim ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
      fprintf( otp_pp , "%lf " , sts_buf_pp[i][j][k] ) ;

    for ( j=0 ; j<Dim  ; j++ )
      for ( k=j+1 ; k<(Dim) ; k++ )
        fprintf( otp_ng , "%lf " , sts_buf_ng[i][j][k] ) ;


    fprintf( otp , "\n" ) ;	
    fprintf( otp_pp , "\n" ) ;
    fprintf( otp_ng , "\n" ) ;
  }

  fclose( otp ) ;
  fclose( otp_pp ) ;
  fclose( otp_ng ) ;

  buff_ind = 0 ;

}*/


extern int t;
FILE *fopen_mkdir(char*, char*);

void write_kspace_data( const char *lbl , complex<double> *kdt ) {
  int i, j , nn[Dim] ;
  FILE *otp ;
  double kv[Dim], k2 ;
  
  char nm[30] ;
  // sprintf( nm, "%s.p.dat" , lbl ) ;

  sprintf( nm, "data/%s.%d.dat" , lbl,t ) ;
  char *mm;
  mm = (char*) nm;

  // char *mode = new char[100];
  // sprintf(mode,"w");
  char *mode= strdup("w");
  char *hehe;
  hehe = (char*) "w";
  otp = fopen_mkdir(mm,hehe);
  // fopen()
  // if(t == 0)
  //   otp = fopen( nm , "w" ) ;
  // else
  //   otp = fopen( nm , "a" ) ;



  // fprintf( otp ,"step %d ,post_spin_dt= %lf\n", t, (step-sample_wait)*delt );

  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
  for ( i=1 ; i<M ; i++ ) {
    unstack( i , nn ) ;

  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
    k2 = get_k( i , kv ) ;
  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
    // cout<< i<<"\t"<<kv<<'\t'<<k2<<endl;
    // cout<<i<<endl;
    // cout<<k2<<endl;

    // for ( j=0 ; j<Dim ; j++ ) 
    //   fprintf( otp , "%lf " , kv[j] ) ;

    fprintf( otp , "%d %1.5e %1.5e %1.5e %1.5e\n" , t, abs(kdt[i]), sqrt(k2), 
       real(kdt[i]) , imag(kdt[i]) ) ;
    // fprintf( otp , "%1.5e %1.5e %1.5e\n" ,abs(kdt[i]), (k2), 
    //     real(kdt[i]) ) ;
    //if ( Dim == 2 && nn[0] == Nx[0]-1 )
      //fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;


}

void write_grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
    
    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}


void rek_mkdir(char *path)
{
  char *sep = strrchr(path, '/' );
  if(sep != NULL) {
    *sep = 0;
    rek_mkdir(path);
    *sep = '/';
  }
  if( mkdir(path,0755) && errno != EEXIST )
    printf("error while trying to create '%s'\n%m\n",path ); 
}


FILE *fopen_mkdir( char *path, char *mode )
{
    char *sep = strrchr(path, '/' );
    if(sep ) { 
       char *path0 = strdup(path);
       path0[ sep - path ] = 0;
       rek_mkdir(path0);
       free(path0);
    } 
    return fopen(path,mode);
}
