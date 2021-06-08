#include "globals.h"
void charge_grid(void);
void write_grid_data(const char*, double*);
void random_config(void);
void write_gro(void);
void read_input_conf(FILE*);

void send_all_x_to_root(void);
int mpi_Bcast_posits(int, int, double**, int, MPI_Comm);
void pa();
void normalize_vector(double * );
void random_vector(double*, bool );
void random_vector(double*);

void initialize_configuration() {
    int i, j;

    if (myrank == 0) {
        random_config();

        FILE* inp;
        inp = fopen("input.lammpstrj", "r");
        if (inp != NULL) {
            read_input_conf(inp);
            fclose(inp);
            printf("input.lammpstrj read!\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (nprocs > 1) {
        mpi_Bcast_posits(nstot, Dim, x, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(mol_type[0]), nstot, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void read_input_conf(FILE* ip) {
    int i, j, k, di;
    double df;
    int rtflag;
    char tt[80];
    fgets(tt, 80, ip);
    fgets(tt, 80, ip);
    fgets(tt, 80, ip);
    fscanf(ip, "%d\n", &di);

    if (di != nstot){
        die("Number of sites in input.lammpstrj does not match!\n");
    }
    printf("\nUsing positions from input.lammpstrj!\n\n");

    fgets(tt, 80, ip);

    fgets(tt, 80, ip);
    fgets(tt, 80, ip);
    fgets(tt, 80, ip);

    fgets(tt, 80, ip);

    for (i = 0; i < nstot; i++) {
        fscanf(ip, "%d", &di);
        fscanf(ip, "%d", &di);
        fscanf(ip, "%d", &di);
        for (j = 0; j < Dim; j++) fscanf(ip, "%lf", &x[i][j]);

        // Read the z position if a 2D simulation
        for (j = Dim; j < 3; j++) fscanf(ip, "%lf", &df);

        fgets(tt, 80, ip);
    }
}

void write_lammps_traj() {
    send_all_x_to_root();

    if (myrank == 0) {
        FILE* otp;
        int i, j, ind, k, m, n, resind;
        if (step == 0){
            otp = fopen("traj.lammpstrj", "w");
        } else {
            otp = fopen("traj.lammpstrj", "a");
        }

        fprintf(otp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n", step,
                nstot);
        fprintf(otp, "ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(otp, "0.0 %lf\n0.0 %lf\n0.0 %lf\n", L[0], L[1],
                (Dim == 3 ? L[2] : 1.0));

        fprintf(otp, "ITEM: ATOMS id type mol x y");
        (Dim >2 ? fprintf(otp," z\n") : fprintf(otp,"\n") );
        ind = 0;
        resind = 0;
        for (k = 0; k < total_num_chains; k++) {
            for (m = 0; m < total_len; m++) {
                fprintf(otp, "%d %d %d", ind + 1, mol_type[ind], resind + 1);

                for (j = 0; j < Dim; j++) fprintf(otp, " %lf", x[ind][j]);

                fprintf(otp, "\n");

                ind++;
            }
            resind++;
        }

        // Outputting the LC only chaines
        ind = nsLCP;
        for (k = 0; k < num_lc_only_c; k++) {
            for (m = 0; m < sites_per_LC; m++) {
                fprintf(otp, "%d %d %d", ind + 1, mol_type[ind], resind + 1);

                for (j = 0; j < Dim; j++) fprintf(otp, " %lf", x[ind][j]);

                fprintf(otp, "\n");

                ind++;
            }
            resind++;
        }
        fclose(otp);

    }  // if ( myrank == 0 )

    MPI_Barrier(MPI_COMM_WORLD);
}

void random_config(void) {
    int i, m, j, k, ind = 0, n, ind_old;
    double rand_u[Dim];
    for (k = 0; k < total_num_chains; k++) {
        // Place the first site of the LC
        ind = k * total_len;
        for (j = 0; j < Dim; j++) {
            x[ind][j] = ran2() * L[j];
        }
        if (Nda > 0) {
            mol_type[ind] = 0;
        } else {
            mol_type[ind] = 1;
        }

        // cout << ind << endl;
            // If the molecules are attached in the middle,
        // ind++;

        // Generate random orientation vector
        random_vector(rand_u, false);

        // Place remaining Polymer A sites along the line
        for (m = 1; m < Nda; m++) {
            ind = k*total_len + m;
            for (j = 0; j < Dim; j++) {
                x[ind][j] = x[ind - 1][j] + ( poly_bond_req + 0.5) * rand_u[j];

                if (x[ind][j] > L[j])
                    x[ind][j] -= L[j];
                else if (x[ind][j] < 0.0)
                    x[ind][j] += L[j];
            }

            mol_type[ind] = 0;

        }
          int cont_ind;

        // Places the polymer types of B down
        if (Nda>0){
           cont_ind = Nda;}
           else {
            cont_ind = 1;
           } 

        for (m = cont_ind; m < polylen; m++) {
            ind = k*total_len + m;
            for (j = 0; j < Dim; j++) {
                x[ind][j] = x[ind - 1][j] + ( poly_bond_req +0.5 ) * rand_u[j];

                if (x[ind][j] > L[j])
                    x[ind][j] -= L[j];
                else if (x[ind][j] < 0.0)
                    x[ind][j] += L[j];
            }

            mol_type[ind] = 1;

        }

        if (side_len <=0) continue;

        // Initializing the side chain

        for (m = 0; m < amount_of_side_chains; m++) {
            ind_old = (k)*total_len + (m)*freq_attached_group;

            random_vector(rand_u, false);

            // If spacers have a length greater than 0
            if (spacer_len > 0) {
                ind = (k)*total_len + polylen + (m)*side_len;

                // Binds the first point of the spacers to the polymer 
            for (j = 0; j < Dim; j++) {
                x[ind][j] = x[ind_old][j] + ( poly_spacer_bond_req +0.1 ) * rand_u[j];
                        if (x[ind][j] > L[j])
                            x[ind][j] -= L[j];
                        else if (x[ind][j] < 0.0)
                            x[ind][j] += L[j];
            }
                // cout << ind << endl;
                    mol_type[ind] = 0;

                //Builds up the spacer chain
                for (n = 1; n < spacer_len ; n++) {
                    ind = (k)*total_len + polylen + (m)*side_len + n ;

                    // cout << ind << endl;
                    for (j = 0; j < Dim; j++) {
                        x[ind][j] =
                            x[ind - 1][j] + ( spacer_bond_req +0.1 ) * rand_u[j];
                        if (x[ind][j] > L[j])
                            x[ind][j] -= L[j];
                        else if (x[ind][j] < 0.0)
                            x[ind][j] += L[j];
                    }
                    mol_type[ind] = 0;
                }
                // ind_old = ind + 1;
                ind_old = ind ;
                mol_type[ind_old] = 0;
            }  // end of spacer section

            if (sites_per_LC <=0 || freq_attached_group ==0) continue;

            ind = (k)*total_len + polylen + (m)*side_len + spacer_len;
            
            // If the molecules are attached in the middle,
            // then update the index. 
            if (not orient_head) {ind =ind + (sites_per_LC-1)/2;}
            
            // If Spacers exist than have the liquid crystaline group bind 
            // to poly
            if (spacer_len > 0) {
                external_lc_bond_req = ( spacer_lc_bond_req +0.1 );
            } else {
                external_lc_bond_req = ( poly_lc_bond_req +0.1 );
            }

            
            for (j = 0; j < Dim; j++) {
              x[ind][j] = x[ind_old][j] + external_lc_bond_req  * rand_u[j];
                if (x[ind][j] > L[j])
                    x[ind][j] -= L[j];
                else if (x[ind][j] < 0.0)
                    x[ind][j] += L[j];
           
            }

            mol_type[ind] = 4;
            random_vector(rand_u);

            if (orient_head) {
                for (n = 1; n < sites_per_LC; n++) {
                    ind = (k)*total_len + polylen + (m)*side_len + spacer_len + n;

                    for (j = 0; j < Dim; j++) {
                        x[ind][j] = x[ind - 1][j] + lc_bond_req * rand_u[j];
                        if (x[ind][j] > L[j])
                            x[ind][j] -= L[j];
                        else if (x[ind][j] < 0.0)
                            x[ind][j] += L[j];
                    }

                    mol_type[ind] = 4;
                }
            } else {

                for (n = 1; n <= (sites_per_LC - 1) / 2; n++) {

                    mol_type[ind] = 4;

                    for (j = 0; j < Dim; j++) {

                        x[ind + n][j] = x[ind + (n - 1)][j] + lc_bond_req * rand_u[j];

                        if (x[ind + n][j] > L[j])
                            x[ind + n][j] -= L[j];
                        else if (x[ind + n][j] < 0.0)
                            x[ind + n][j] += L[j];

                        x[ind - n][j] = x[ind - (n - 1)][j] - lc_bond_req * rand_u[j];

                        if (x[ind - n][j] > L[j])
                            x[ind - n][j] -= L[j];
                        else if (x[ind - n][j] < 0.0)
                            x[ind - n][j] += L[j];
                    }
                    mol_type[ind + n] = 4;
                    mol_type[ind - n] = 4;
                }

            } // end of lc_section
        }
    }  // for n=0:Nd-1



    // PURE LC SECTION
      ind = nsLCP;
  for ( k=0 ; k< num_lc_only_c ; k++ ) { 

    // Place the first site of the LC
    for ( j=0 ; j<Dim ; j++ ) x[ind][j] = ran2() * L[j] ;
    mol_type[ind] = 4;

    ind++ ;
      
    // Generate random orientation vector
    double rand_u[Dim];
    random_vector(rand_u);


    // Place remaining sites in a line along u
    for ( m=1 ; m<sites_per_LC; m++ ) {
    
      for ( j=0 ; j<Dim ; j++ ) {
        x[ind][j] = x[ ind-1 ][ j ] + ( lc_bond_req +0.1 ) * rand_u[j];
 
        if ( x[ind][j] > L[j] )
          x[ind][j] -= L[j] ;
        else if ( x[ind][j] < 0.0 )
          x[ind][j] += L[j] ;
      }
 
      mol_type[ind] = 4;
      ind++ ;
      
    }
  }// for n=0:Nd-1

    printf("Random config generated!\n");
    fflush(stdout);
}



void pa(){
    for (int i = 0; i < nstot; i++){
        cout << "particle " << i << ": "; 
        for (int j = 0; j<Dim; j++){
           cout << x[i][j] << " " ;
        }
        cout << endl;
    }
}

void random_vector(double* vector, bool force_dir){
    if (force_dir){
        for (int j = 0; j < Dim; j++) vector[j] = 1.0;
        normalize_vector(vector);
    }
    else {
        for (int j = 0; j < Dim; j++) vector[j] = gasdev2();
        normalize_vector(vector);   
    }
}

void random_vector(double* vector){
    random_vector(vector, init_flag == 1);
}
