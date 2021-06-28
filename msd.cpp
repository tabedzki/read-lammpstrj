#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include "log_space.h"

using namespace std;

void log_space(int,int,int);
double pbc_mdr2(vector<double>, vector<double>, vector<double>, vector<double>);

void calc_msd( vector<vector<vector<double>>> xt,
    const int nsites, 
    const int nfr, 
    vector<vector<double>> L ) {

  log_space(1,nfr,nfr/2) ;

  vector<double> msd(n_times,0), dr(3);

  for ( int ti = 0; ti < n_times ; ti++ ) {

    int delt = times[ti] ;

    for ( int t=0 ; t<nfr-delt ; t++ ) {

      for ( int i=0 ; i<nsites ; i++ ) {

        double mdr2 = 0.0;
        for ( int j=0 ; j<3 ; j++ ) {
          dr[j] = xt[t+delt][i][j] - xt[t][i][j] ;
          mdr2 += dr[j] * dr[j];
        }

        msd[ti] += mdr2 ;

      }// i=0:nsites

    }// t=0:nfr-delt

  }// ti=0:n_times

  ofstream otp("msd.dat");

  for ( int ti=0 ; ti<n_times ; ti++ ) {
    int delt = times[ti] ;

    msd[ti] *= ( 1.0 / double(nsites * (nfr-delt) ) ) ;

    otp << delt << " " << msd[ti] << endl;

  }

  otp.close() ;
}