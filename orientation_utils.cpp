#include "globals.h"

void normalize_vector(double * );
void rg_order();
void pbc_mdr3(double*, double*, double*);
void shift_poly_BB_into_one_box(int , int , double*);
void shift_poly_BB_into_one_box(int , int , double*[Dim]);
void shift_poly_BB_into_one_box(Array<double, -1 ,Dim> & );
void shift_poly_BB_into_one_box(Map<Array<double, -1 ,Dim>> & );
Matrix<double,Dim,Dim> Eye = Matrix<double,Dim,Dim>::Identity();
extern double Rg2;

// Takes the particle-level S tensors and adds them to the grid //
// It assumes that particle_orientations has been calculated    //
// prior to calling this function.
void add_S_to_grid( ) { 

  if (sites_per_LC == 0 ) return;

  int i,j,k,m, id, Mid, mod ;
  double W3 ;

  
  for ( i=0 ; i<ns_loc+total_ghost ; i++ ) {
    if ( i < ns_loc ) 
      id = my_inds[i] ;
    else
      id = ghost_inds[ i-ns_loc ] ;

    if ( !LC_has_orientation( id ) )
      continue ;


    for ( m=0 ; m<grid_per_partic ; m++ ) {
      Mid = grid_inds[id][m] ;

      if ( Mid == -1 )
        continue ;

      W3 = grid_W[id][m] ;

      for ( j=0 ; j<Dim ; j++ )
        for ( k=0 ; k<Dim ; k++ ){ 
          S_field[Mid][j][k] += W3 * mono_S[id][j][k] ;
          }
    }// for m < grid_per_partic

  }// for i < ns_loc + total_ghost

}



// Calculates the u vector 
// for each LC molecule
void particle_orientations(bool local) {

  if (sites_per_LC == 0 ) return;

  int i,j,id,k,m,n, mod,mod1, id1 ;
  double mdr, mdr2, dr[Dim], drij[Dim] ;
  int total = 0, l_total= 0;
  int values; 
  (local ? values = ns_loc : values = nstot );
  for ( i=0 ; i< values; i++ ) {
    // id = my_inds[i] ;
    (local ? id = my_inds[i] : id = i);

    if ( !LC_has_orientation(id) ) continue ;
    l_total++;
    // Reset monomer orientation
    for ( m=0 ; m<Dim ; m++ ){
      mono_u[id][m] = 0.0 ;
    }

    // Orientation vector with the next monomer
    id1 = id + 1 ;
 
    // Define unit vector
    mdr2 = pbc_mdr2( x[id], x[id1], dr ) ;
    mdr = sqrt( mdr2 ) ;
    for ( m=0 ; m<Dim ; m++ ) dr[m] /= mdr ;

    // Accumulate u, S
    for ( m=0 ; m<Dim ; m++ ) mono_u[id][m]  = dr[m] ;
 
    // Orientation vector with the previous monomer
    id1 = id - 1 ;
 
    // Define unit vector
    mdr2 = pbc_mdr2( x[id1], x[id], drij ) ;
    mdr = sqrt( mdr2 ) ;
    for ( m=0 ; m<Dim ; m++ ) mono_u[id][m]  += drij[m]/mdr;
    normalize_vector(mono_u[id]);

  }// for i < ns_loc

}

void particle_orientations(int j){
  if (j==1){
    particle_orientations(true);
  }
  else{
    particle_orientations(false);
  }
}

void particle_orientations(void ){
  particle_orientations(true);
}

// Calculates the S tensor for each particle.
// Assumes u has been calculated and communicated
// in parallel across adjacent procs.
void particle_Stensor() {
  if (sites_per_LC == 0 ) return;

    int i, j, id, k, m, n, mod, mod1, id1;
    double mdr, mdr2, dr[Dim];

    for (i = 0; i < ns_loc; i++) {
        id = my_inds[i];

        if (!LC_has_orientation(id)) continue;

        for (m = 0; m < Dim; m++) {
            for (n = 0; n < Dim; n++) {
                mono_S[id][m][n] =
                    mono_u[id][m] * mono_u[id][n] - KDelta(m, n) / Dim;
            }
        }
    }  // i=0:ns_loc
}

void gubbins_q_tensor(){
  
  if (sites_per_LC == 0 ) return;

    int i, j, id, k, m, n ;
    double mdr, mdr2, dr[Dim];
    double  total=0;
    double l_total=0;
    Matrix<double,Dim,Dim> Q_ten = Matrix<double,Dim,Dim>::Zero() ;
    Matrix<double,Dim,Dim> Q_ten_loc = Matrix<double,Dim,Dim>::Zero() ;
    EigenSolver<Matrix<double,Dim,Dim>> ces ;
    Array<double,Dim,1> Eigen_vals;
    Matrix<double,Dim,1> mono_u_loc;

    for (i = 0; i < ns_loc; i++) {
        id = my_inds[i];

        if (!LC_has_orientation(id)) continue;
        l_total++;

        for (m=0;m<Dim;m++) mono_u_loc(m)=mono_u[id][m];
        Q_ten_loc += mono_u_loc * mono_u_loc.transpose();

    }  // i=0:ns_loc


    MPI_Allreduce(Q_ten_loc.data(), Q_ten.data(), Dim*Dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD  );
    MPI_Allreduce( &l_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

    Q_ten /= total;
    Q_ten =1.0 /(static_cast<double>(Dim) - 1.0 ) * 
        (static_cast<double>( Dim )*Q_ten - Eye);
    ces.compute(Q_ten);

    Eigen_vals= ces.eigenvalues().array().real().cast<double>();
    gubbins_max_real = Eigen_vals.maxCoeff(&m); 
    Eigen_vals.abs().maxCoeff(&m); 
    gubbins_largest= Eigen_vals(m);
    director = ces.eigenvectors().col(m).real().cast<double>() ;
}

void calculate_polymer_segments(void ){
  int index;
  if (myrank==0){
  for (int i=0; i< total_num_chains; i++ ){
    for (int j=0; j<(polylen-1); j++){
      index = i*total_len + j;  
      pbc_mdr2(x[index+1], x[index], dr_store[index] ); // change this to stored
    }
  }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}




/**
 * @brief Only call this after sending all x to root
 * 
 */
void rg_order(){
  int index, m;
  if (myrank == 0){
    double eigen_val;
    Matrix<double,Dim,Dim> Rg_tensor = Matrix<double,Dim,Dim>::Zero();
    EigenSolver<Matrix<double,Dim,Dim>> ces ;
    Matrix<double, Dim, 1> local_director = Matrix<double,Dim,1>::Zero();
    Matrix<double,1,Dim> temp;
    poly_lc_order = 0;
    Rg2 = 0;
    Matrix <double, -1, 1> nm_dotted;


    Array<double, -1 ,Dim> local_backbone = Array<double, -1 ,Dim>::Zero(polylen,Dim);

    for(int i=0; i < total_num_chains; i++){ 
      index= i * total_len;

      for (int j=0; j<polylen; j++){
        for (int dim=0; dim < Dim; dim ++) local_backbone(j,dim) = x[j+index][dim];
      }
      shift_poly_BB_into_one_box(local_backbone);
      local_backbone = local_backbone.rowwise() - local_backbone.colwise().mean();

      Rg_tensor.Zero();  
      for (int j = 0; j < polylen ; j++) {
        temp = local_backbone.row(j).matrix();
        Rg_tensor += temp.transpose() * temp;
        Rg2 += temp.dot(temp);
      }

      Rg_tensor = 1.0/ (Dim - 1.0) * 
        (static_cast<double>( Dim ) * Rg_tensor/polylen - Eye);

      ces.compute(Rg_tensor);
      eigen_val = ces.eigenvalues().real().maxCoeff(&m);
      local_director = ces.eigenvectors().real().col(m);
      double dot_value = (local_director.dot(director));
      poly_lc_order += 
        1.0/(2.0) * 
        ( 3.0  * dot_value*dot_value - 1);


    }
    Rg2 /= polylen;
    Rg2 /= total_num_chains;
    poly_lc_order /= total_num_chains;

  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void doi_q_tensor(){
  
  if (sites_per_LC == 0 ) return;
  
    int i,  id,  m;
    double mdr,  summed_dot=0;
    int total=0;
    int l_total=0;
    double ndir_gather;
    double n_dir_local[Dim];
    for (m=0; m<Dim;m++) n_dir_local[m]=0;

    for (i = 0; i < ns_loc; i++) {
        id = my_inds[i];

        if (!LC_has_orientation(id)) continue;
        l_total++;
        for (m = 0; m < Dim; m++) {
          n_dir_local[m] += mono_u[id][m];
        }

    }  // i=0:ns_loc

    MPI_Allreduce(n_dir_local, n_dir, Dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&l_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (m = 0; m < Dim; m++) {
      n_dir[m] /= total;
    }
    normalize_vector(n_dir);
    l_total=0;

    for (m=0; m<Dim;m++) n_dir_local[m]=0;

    for (i = 0; i < ns_loc; i++) {
      id = my_inds[i];

      if (!LC_has_orientation(id)) continue;
      l_total++;
      mdr = 0; 
      for( m = 0; m < Dim ; m++) mdr+= mono_u[id][m] * n_dir[m]; 
       summed_dot += mdr * mdr;
    }

    // MPI_Allreduce(&ndir_tmp, &ndir, Dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD  );
    MPI_Allreduce( &l_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

    MPI_Allreduce( &summed_dot, &ndir_gather, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
    doi_lc= Dim/(Dim-1.0)*  (ndir_gather/(total ) - 1./Dim);
    
}

// Determines if a particle carries a u vector and
// and S tensor.
int LC_has_orientation( int id ) {
  if (mol_type[id] != 4) return 0;
  if ((mol_type[id-1] != 4) || mol_type[id+1] != 4  ) return 0;
  int id_old =id;
  if (total_len !=0){
  id %= total_len;
  id -= polylen;
  id %= side_len;
  id -= spacer_len;

  // if ( ( id%sites_per_LC ) == ( sites_per_LC/2 ) )
  if (id == (sites_per_LC-1)/2 ){
    // cout << id_old << endl;
    return 1;
    }
  }
  if (not(nstot == nsLCP)) {
  id = id_old;
  id -= total_num_chains * total_len;
  id %= sites_per_LC; 
    if (id == (sites_per_LC - 1) / 2) {
    // cout << id_old << endl;
    return 1;
  }
  }
  return 0;
}


void make_eigenval_map() {

  int i,j ;
  double *Vec ;
  Vec = new double [Dim] ;

  // cout << 1 << endl;
  for ( i=0 ; i<ML ; i++ ) {
    tmp[i] = diag_mat( S_field[i], Vec ) ;
  }

  // cout << 2 << endl;
  delete Vec ;

  lc_order_param = integrate(tmp) / V ;

  write_grid_data( "eigen_val", tmp ) ;

}

// Calculates the convolution of the S field with the Gaussian potential.
// This is needed for both force and energy evaluations. 
// Should be calculated at the beginning of the force routine.
void calc_S_conv_uG() {
    int i, j, k, m, n;

    fftw_fwd(uG, ktmp2);

    for (m = 0; m < Dim; m++) {
        for (n = 0; n < Dim; n++) {
            for (i = 0; i < ML; i++) tmp[i] = S_field[i][m][n];

            fftw_fwd(tmp, ktmp);

            for (i = 0; i < ML; i++) ktmp[i] *= ktmp2[i];

            fftw_back(ktmp, tmp);

            for (i = 0; i < ML; i++) S_conv_u[i][m][n] = tmp[i];
        }
    }
}

void normalize_vector(double* t1 ){
  double sum=0;
  int m;
  for (m=0;m<Dim; m++) sum+=t1[m]*t1[m];
  sum = sqrt(sum);
  for (m=0;m<Dim; m++){t1[m]/= sum;}

}

void zero_vector(double* t1 ){
  int m;
  for (m=0;m<Dim; m++) t1[m] =0;

}

double P2(double x){
  return (1.5*x*x - 0.5);
}
