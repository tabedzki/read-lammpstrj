#include "globals.h"
void bond_stress();
void nonbond_stress();

void calc_stress() {


  int j1, j2 ;
  nonbond_stress() ;
  bond_stress();

  for ( j1=0 ; j1<Dim ; j1++ )
    for ( j2=0 ; j2<Dim ; j2++ )
      //Ptens[j1][j2] = -(Stress_nb[j1][j2]) + Stress_bonds[j1][j2] + rho0 * Kdelta(j1,j2) ;
      //Ptens[j1][j2] = -(Stress_nb[j1][j2] + Stress_ss[j1][j2]) + Stress_bonds[j1][j2]  + rho0 * Kdelta(j1,j2) ;
      // Ptens[j1][j2] = -(Stress_nb[j1][j2]) ;
      Ptens[j1][j2] = -(Stress_nb[j1][j2] + Stress_bonds[j1][j2]) + nD/V * KDelta(j1,j2);
      // Ptens[j1][j2] = -(Stress_nb[j1][j2] ) ;
      //Ptens[j1][j2] = -Stress_nb[j1][j2]
        //+ double(nA)/V * Kdelta(j1,j2) ;

  Pscalar = 0.0 ;
  for ( j1=0 ; j1<Dim ; j1++ )
    // Pscalar += Rg3 * Ptens[j1][j1] / double( Dim ) + Rg3*nstot/V /3;
    // Pscalar += Rg3 * Ptens[j1][j1] / double( Dim ) + Rg3*nstot/V/3 ;
    Pscalar += Ptens[j1][j1] / double( Dim );

}


void nonbond_stress() {

  int i, j1, j2 ;

      // for ( i=0 ; i<M ; i++ )
      //   rhot[i] -= rho0;
  fftw_fwd( rhot , ktmp ) ;
   for(i=0; i<ntypes; i++)
	fftw_fwd( rho[i] , rho_hat[i] ) ;
  // fftw_fwd( rhoga,tmp_Ng);

 

  for ( j1=0 ; j1<Dim ; j1++ ) {
    for ( j2=0 ; j2<Dim ; j2++ ) {
    //compressibility 
     if(kappa>0){
      for ( i=0 ; i<M ; i++ ){
        ktmp2[i] =  ktmp[i] * vir_func_hat[j1][j2][i] ;
      }
      fftw_back( ktmp2 , tmp ) ;

      for ( i=0 ; i<M ; i++ ){
        tmp2[i] = (rhot[i] )* tmp[i]*kappa/2.0/V ; 
  	  }

    }


    //  // Flory-Huggins
    //  if(chiAB>0){
    //   for ( i=0 ; i<M ; i++ ){
    //     ktmp2[i] = ( rho_hat[1][i]) * vir_func_hat[j1][j2][i] ;
        
    //   }

    //   fftw_back( ktmp2 , tmp ) ;
    
    //   for ( i=0 ; i<M ; i++ )
    //     tmp2[i]  += (rho[0][i] )* tmp[i] *chiAB/V/rho0; 
    //  }

      //particle-particle contribution 
    //  if(nP >0 ){      
    //   if(nP>1){
    //     for ( i=0 ; i<M ; i++ ){
          
    //     ktmp2[i] = ( rho_hat[2][i]) * vir_funcpp_hat[j1][j2][i] ;
            
    //     }

    //     fftw_back( ktmp2 , tmp_PP ) ;
    //     for ( i=0 ; i<M ; i++ ){
    //       tmp_PP[i]  *= (rho[2][i])* kappa_p/V/2.0/rho0; 
    //       tmp2[i] += tmp_PP[i];

    //     }
    //     Stress_PP[j1][j2] = integrate( tmp_PP )  ;

    //   }//nP >1

    //   // partilce-polymer contribution
    //   for ( i=0 ; i<M ; i++ ){
    //     ktmp2[i] = ( rho_hat[2][i]) * vir_funcpg_hat[j1][j2][i] ;
    //   }

    //   fftw_back( ktmp2 , tmp ) ;
    //   for ( i=0 ; i<M ; i++ )
    //     tmp2[i]  += (rho[0][i]*kappa + (A_partics > 0? chiAB+kappa: kappa)*rho[1][i])* tmp[i] /V/rho0; 

    //  }//nP >0 

      Stress_nb[j1][j2] = -integrate( tmp2 )  ;
    }
  }

}


void bond_stress() {
  int center_ind,i,k,j1, j2, m, ind ;
  double mdr,mdr2, dr[Dim] ;

  ind = 0 ;



  for (j1=0 ; j1<Dim ; j1++ ) 
    for ( j2=0 ; j2<Dim ; j2++ )
      Stress_bonds[j1][j2] = 0.0 ;
	// for(m=0; m<nthreads; m++)
	// 	 Stress_bond_t[j1][j2][m] = 0.0 ;
  //   


  // Diblock bonds //
// #pragma omp parallel for private(ind,m,mdr2,dr,j1,j2)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      // int tid = omp_get_thread_num() ;
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
          Stress_bonds[j1][j2] +=  dr[j1] * dr[j2] ;
      
    }
  
  } // for ( i=0 ; i<nT[k]


   // Homopolymer A bonds //
 // #pragma omp parallel for private(ind,m,dr,mdr2,j1,j2)
   /* for ( i=0 ; i<nA ; i++ ) { */
   /*   for ( m=0 ; m<Nha - 1 ; m++ ) { */

   /*     ind = nD * (Nda + Ndb) + i * Nha + m ; */

   /*     mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ; */
   /*     int tid = omp_get_thread_num() ; */

   /*     for ( j1=0 ; j1<Dim ; j1++ ) */
   /*       for ( j2=0 ; j2<Dim ; j2++ ) */
   /*         Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ; */
    
   /*   } */
  
//   } // for ( i=0 ; i<nT[k]

//   // Homopolymer B bonds //
// // #pragma omp parallel for private(ind,m,dr,mdr2,j1,j2)
//   for ( i=0 ; i<nB ; i++ ) {
//     for ( m=0 ; m<Nhb - 1 ; m++ ) {

//       ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

//       mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;
   
//       int tid = omp_get_thread_num() ;

//       for ( j1=0 ; j1<Dim ; j1++ )
//         for ( j2=0 ; j2<Dim ; j2++ )
//           Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;
      
//     }

//   } // for ( i=0 ; i<nT[k]

// // HomoA gaft chain bons //
//   // #pragma omp parallel for \
//   // private( dr,center_ind, m, ind, mdr2, j1,j2, k, prev_graft_site ) 
//   for ( i=0 ; i<nP ; i++ ) {

//     center_ind = nD * ( Nda + Ndb ) + nA * Nha + nB * Nhb + sites_per_gnp * i ;

//     for ( m=0 ; m<ng_per_partic ; m++ ) {
//       ind = center_ind + m * ( Ng + 1 ) + 1 ;

//       int tid = omp_get_thread_num() ;
 
//       prev_graft_site = ind ;
//       ind += 1 ;
//       for ( k=0 ; k<Ng ; k++ ) {
//         mdr2 = pbc_mdr2( x[ind-1] , x[ind] , dr ) ;

// 	for ( j1=0 ; j1<Dim ; j1++ )
//         	for ( j2=0 ; j2<Dim ; j2++ )
//           		Stress_bond_t[j1][j2][tid] +=  dr[j1] * dr[j2] ;

//     	  ind++;
// 	}

//       }//m

//   }//for ( i=0 ; i<nP 


  // for ( j1=0 ; j1<Dim ; j1++ )
	// for ( j2=0 ; j2<Dim ; j2++ )
	//  for ( m=0 ; m<nthreads ; m++ )
	// 	Stress_bonds[j1][j2] += 3.0 * Stress_bond_t[j1][j2][m] /V;

}
