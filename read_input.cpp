#include "globals.h"
#include "stiff-polymer-backbone.h"


void read_anneal( void ) ;

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[180] ;

  inp = fopen( "lc.input" , "r" ) ;
  if ( inp == NULL ) 
    die("Failed to open lc.input!");
    
  // Skip header
  fgets(tt, 80, inp) ;
  
  // Lengths of polymer
  // fscanf( inp, "%d", &sites_per_LC ) ;
  fscanf( inp, "%d", &Nda ) ;
  fscanf( inp, "%d", &spacer_len ) ;
  fscanf( inp, "%d", &sites_per_LC ) ;
  fscanf( inp, "%d", &Ndb ) ;

  SHOW(Nda);
  SHOW(spacer_len);
  SHOW(sites_per_LC);
  SHOW(Ndb);
  sites_per_LC = sites_per_LC ; //Restoring the proper value 

  fgets(tt, 80, inp) ;

  //Freq
  fscanf( inp , "%d" , &freq_attached_group) ; 
  fgets(tt, 80, inp) ;
  SHOW(freq_attached_group);

  fscanf( inp , "%lf" , &a_squared  ) ; 
  SHOW(a_squared);
  a_squared *= a_squared ;
  fgets( tt , 80 , inp ) ;

  fscanf(inp, "%lf %lf", &poly_bond_req, &spacer_bond_req ) ;
  fscanf(inp, "%lf", &lc_bond_req  ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(poly_bond_req);
  SHOW(spacer_bond_req);
  SHOW(lc_bond_req);

  fscanf(inp, "%lf %lf", &poly_bond_k, &spacer_bond_k ) ;
  fscanf(inp, "%lf", &lc_bond_k  ) ;
  SHOW(poly_bond_k);
  SHOW(spacer_bond_k);
  SHOW(lc_bond_k);
  fgets( tt , 80 , inp ) ;

  fscanf(inp, "%lf %lf", &poly_spacer_bond_req, &spacer_lc_bond_req ) ;
  fscanf(inp, "%lf", &poly_lc_bond_req ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(poly_spacer_bond_req);
  SHOW(spacer_lc_bond_req);
  SHOW(poly_lc_bond_req);

  fscanf(inp, "%lf %lf", &poly_spacer_bond_k, &spacer_lc_bond_k ) ;
  fscanf(inp, "%lf", &poly_lc_bond_k ) ;
  SHOW(poly_spacer_bond_k);
  SHOW(spacer_lc_bond_k);
  SHOW(poly_lc_bond_k);
  fgets( tt , 80 , inp ) ;

  // Harmonic Angle Force Constance
  fscanf( inp , "%lf", &lc_angle_k ) ;
  SHOW(lc_angle_k);
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;

  //////////////////////////////////////////
  // Homopolymer based parameters and chi //
  //////////////////////////////////////////
  fscanf( inp , "%lf %d" , &phiHA , &Nha ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(phiHA); SHOW(Nha);

  fscanf( inp , "%lf %d" , &phiHB , &Nhb ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(phiHB); SHOW(Nhb);

  fscanf( inp , "%lf %d" , &phiHC , &Nhc ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(phiHC); SHOW(Nhc);

  fscanf( inp , "%lf %lf %lf" , &chiAB, &chiAC, &chiBC ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(chiAB); SHOW(chiAC); SHOW(chiBC);

  ////////////////////////////////////////
  // Shearable bendable model parameters//
  ////////////////////////////////////////

  // Blank line //
  fgets( tt , 80 , inp ) ;
  { using namespace stiff;
  fscanf( inp , "%lf %lf %lf" , &epsilon_b, &epsilon_parallel, &epsilon_perp );
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%lf %lf %lf" , &eta, &stiff::gamma, &l0);
  fgets( tt , 80 , inp ) ;
}

  // Skipping gauss
  fscanf( inp, "%d %d %d %d", &do_gauss, &do_wlc, &do_field, &do_hard_angles_LC);
  do_dsswlc= 0;
  fgets( tt , 80 , inp ) ;

  // Blank line //
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%lf" , &phiP ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(phiP);

  fscanf( inp , "%lf" , &Rp ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(Rp);

  fscanf( inp , "%lf" , &Xi ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(Xi);

  fscanf( inp , "%d" , &A_partics ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(A_partics);

  ///////////////////////////
  // Simulation parameters //
  ///////////////////////////
  fscanf( inp , "%lf" , &rho0 ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(rho0);

  fscanf( inp , "%lf" , &phiFLC) ;
  fgets( tt , 80 , inp ) ;
  SHOW(phiFLC);

  fscanf( inp , "%lf" , &kappa ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(kappa);

  fscanf( inp , "%lf" , &mu ) ;
  SHOW(mu);
  fgets( tt , 80 , inp ) ;
  if (sites_per_LC == 0) mu = 0;

  // Orient the LC group by the head
  fscanf( inp , "%d" , &orient_head ) ;
  SHOW(orient_head);
  fgets( tt , 80 , inp ) ;

  Diff = ( double* ) calloc( 5 , sizeof( double ) ) ;
  fscanf( inp , "%lf %lf %lf" , &Diff[0] , &Diff[1], &Diff[3] ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(Diff[0]); SHOW(Diff[1]); SHOW(Diff[3]);
  if (myrank==0) cout << "Diff of LC is set to be that of polymers."<< endl;
  Diff[4] = Diff[0];
  // Blank line //
  fgets( tt , 80 , inp ) ;

  VectorXf Lengths(Dim);
  VectorXf Lengths_half(Dim);

  // Lengths //
  for ( i=0 ; i<Dim ; i++ ) {
    fscanf( inp , "%lf" , &L[i] ) ; 
    Lengths(i) = L[i];
    }
    if (myrank==0) cout << Lengths << endl;
    Lengths_half = Lengths/2;
  fgets( tt , 80 , inp ) ;
  cout <<"Lengths " << Lengths <<endl;
  // Grid dimensions //
  for (i = 0; i < Dim; i++) {
      fscanf(inp, "%d", &Nx[i]);
      if (Nx[i] % 2 == 0) {
          if (myrank == 0) {
              cout << "WARNING: Your gridpoints have increased from " << Nx[i] << " to " << (Nx[i] + 1);
              cout << " for dimension " << i << "." << endl;
              cout << "Please, in the future, use a number that is odd or a multiple of odd prime numbers." << endl;
          }
          Nx[i]++;
      }
    SHOW(Nx[i]);
  }
  fgets( tt , 80 , inp ) ;


  // Time step//
  fscanf( inp , "%lf", &delt ) ;
  SHOW(delt);
  fgets( tt , 80 , inp ) ;

  // Grid interpolation order //
  fscanf( inp , "%d" , &pmeorder ) ;
  SHOW(pmeorder);
  fgets( tt , 80 , inp ) ;

  // Buffer between parallel domains //
  fscanf( inp , "%lf" , &send_buff ) ;
  SHOW(send_buff);
  fgets( tt , 80 , inp ) ;


  // Blank line //
  fgets( tt , 80 , inp ) ;


  fscanf( inp , "%d" , &nsteps ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d" , &print_freq ) ;
  fgets( tt , 80 , inp ) ;
  fscanf( inp , "%d %d" , &sample_wait , &sample_freq ) ;
  fgets( tt , 80 , inp ) ;
  SHOW(nsteps);
  SHOW(print_freq);
  SHOW(sample_wait);
  SHOW(sample_freq);
  
  fscanf( inp , "%d" , &traj_freq ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d" , &init_flag ) ;
  fgets( tt , 80 , inp ) ;
  
  fscanf( inp , "%d" , &print_vectors) ;
  fgets( tt , 80 , inp ) ;

  fscanf( inp , "%d" , &multi_vec_files) ;
  fgets( tt , 90 , inp ) ;

  fscanf( inp , "%d" , &print_order_params) ;
  fgets( tt , 80 , inp ) ;

  read_anneal() ;

}
