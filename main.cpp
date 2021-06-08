#define MAIN
#include "globals.h"
#include "timing.h"
#include "stiff-polymer-backbone.h"

void update_positions( void ) ;
void send_all_x_to_root(void);
void initialize( void ) ;
void write_lammps_traj( void ) ;
void write_vector_fields( void ) ;
void write_grid( void ) ;
void forces( void ) ;
double integrate( double* ) ;
void write_stress( void ) ;
void calc_Unb( void ) ;
void anneal_update( void ) ;
void make_eigenval_map( void ) ;
void gubbins_q_tensor();
void doi_q_tensor();
void particle_orientations();
void particle_orientations(bool);
void rg_order();
void calculate_polymer_segments();
void calc_scalar_parameters();
double calculate_R_parallel();
double calculate_Rg2();
void print_screen();
void print_log(FILE*);
bool prepare_rg_variables(bool);
void calc_stress();

double R_parallel,R_parallel_sq, R_perpendicular_sq, R_perpendicular, Rg2;
extern double sp_plain;


int main( int argc , char** argv ) {

  int i,j,k ;
  
  main_t_in = time(0) ;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank ) ;
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs ) ;
  fftw_mpi_init() ;

  initialize() ;

  write_lammps_traj();

  prepare_rg_variables(false);
    if (myrank == 0) {
      cout << "Rg2 " << Rg2 << " R||2 " << R_parallel_sq << " Rp2 "
           << R_perpendicular_sq << " aniso " << anisotropic_ratio
           << " sp: " << sp_plain
           << " alt_calc: " << sqrt((1 + 2 * sp_plain) / (1 - sp_plain))
           << endl;
    }

  if (print_vectors == 1) {
    if (myrank == 0) particle_orientations(false);
    write_vector_fields();
  }

  FILE* otp;
  if (myrank == 0) {
    otp = fopen("data.dat", "w");
    fprintf(otp, "step Ubond Uangle Ukappa ");
    fprintf(otp, "Ubend Ucomp Ushear ");
    // if (mu != 0.0 || print_order_params == 1){fprintf(otp , "U_ms
    // lc_order_param doi_param gubbins_max gubbins_largest");};
    if (mu != 0.0 || print_order_params == 1) {
      fprintf(otp, "U_ms lc_order_param gubbins_max gubbins_largest");
      if (not polylen == 0.0) fprintf(otp, " poly_lc_order");
      if (not polylen == 0.0) fprintf(otp, " aniso_ratio");
    }
    fprintf(otp, "\n");
  }

  if (myrank == 0) {
    printf("Entering main loop!\n");
    fflush(stdout);
  }

  calc_Unb();
  if (mu != 0.0 || print_order_params == 1) {
    make_eigenval_map();
    // doi_q_tensor();
    gubbins_q_tensor();
  }

  print_screen();

  fftw_fwd(rho[0], ktmp);
  for (i = 0; i < ML; i++) ktmp[i] = ktmp[i] * conj(ktmp[i]);
  write_kspace_data("sk", ktmp);

  if (mu != 0.0) write_Sfield_data("Sfield", S_field);

  if (nLC > 0.0) write_grid_data("rholc", rholc);

  if (step > sample_wait) {
    for (i = 0; i < ML; i++) ktmp2[i] = avg_sk[0][i] / num_averages;

    write_kspace_data("avg_sk", ktmp2);
      }

      print_log(otp);


      bool rg_done ;
      for (step = 1; step <= nsteps; step++) {

        rg_done = false;

        if (do_anneal && step == next_anneal_update) anneal_update();

        ////////////////////
        // Core algorithm //
        ////////////////////
        forces();
        update_positions();

        ////////////////////////////////
        // Calculate structure factor //
        ////////////////////////////////
        if (step > sample_wait && step % sample_freq == 0) {
          fftw_fwd(rho[0], ktmp);
          for (i = 0; i < ML; i++) {
            double k2, kv[Dim];
            k2 = get_k(i, kv);
            avg_sk[0][i] += ktmp[i] * conj(ktmp[i]);
          }
          num_averages += 1.0;
        }

        /////////////////////////
        // Write lammps output //
        /////////////////////////
        if (step % traj_freq == 0 && traj_freq > 0) {
          write_lammps_traj();
          rg_done = prepare_rg_variables(rg_done);
          if (not polylen == 0.0) {
            if (myrank==0)
                cout << "Rg2 " << Rg2 << " R||2 " << R_parallel_sq<< " Rp2 " << R_perpendicular_sq << " aniso " << anisotropic_ratio <<
                " sp: " << sp_plain  << " alt_calc: " << sqrt((1+2*sp_plain)/(1-sp_plain)) << endl;
          }

          // SHOW(director.transpose());
          if (print_vectors == 1) {
            write_vector_fields();
          }
        }

        ///////////////////
        // Write outputs //
        ///////////////////
        if (step % print_freq == 0 || step == nsteps - 1) {
          calc_Unb();
          if (mu != 0.0 || print_order_params == 1) {
            make_eigenval_map();
            // doi_q_tensor();
            gubbins_q_tensor();
          }
          rg_done = prepare_rg_variables(rg_done);

          print_screen();

          fftw_fwd(rho[0], ktmp);
          for (i = 0; i < ML; i++) ktmp[i] = ktmp[i] * conj(ktmp[i]);
          write_kspace_data("sk", ktmp);

          if (mu != 0.0) write_Sfield_data("Sfield", S_field);
          if (nLC > 0.0) write_grid_data("rholc", rholc);
          if (step > sample_wait) {
            for (i = 0; i < ML; i++) ktmp2[i] = avg_sk[0][i] / num_averages;

            write_kspace_data("avg_sk", ktmp2);
          }
          print_log(otp);

        }  // if step % print_Freq == 0
  }// for step=0:max_steps

  if ( myrank == 0 ) 
    fclose( otp ) ;

  main_t_out = time(0); 
  if ( myrank == 0 ) {
    printf("Total time: %d mins, tot seconds: %ds\n", (main_t_out - main_t_in)/60, 
        (main_t_out-main_t_in) ) ;
 
    printf("FFT time: %d mins %d secs\n", fft_tot_time/60, fft_tot_time%60 ) ;
    printf("Update time: %dmins %d secs\n", move_tot_time/60, move_tot_time%60 ) ;
    printf("Grid time: %dmins %d secs, tot seconds: %ds\n", grid_tot_time/60, grid_tot_time%60, grid_tot_time ) ;

    if ( nprocs > 1 ) {
      printf("\nUpdate time, not including comm: %d mins %d secs\n", move_minus_comm/60, move_minus_comm%60);
      printf("Total time spent in swap/comm routines: %d mins %d secs, tot sec: %d\n", swap_tot_time/60 , swap_tot_time%60, swap_tot_time ) ;
      printf("Total time for exchanging forces: %d s, ghosts: %d s, partics: %d s\n",
          force_comm_tot_time, swap_ghosts_tot_time, swap_partics_tot_time );
      if ( time_debug_tot_time > 0 ) 
        printf("\nDebug total time: %d\n", time_debug_tot_time ) ;
    }
  }

  fftw_mpi_cleanup();

  MPI_Finalize() ;

  return 0 ;

}


void print_screen(){

  if (myrank == 0) {
    printf("step %d of %d  Ubond: %lf ", step, nsteps, Ubond);
    printf("Uangle: %lf ", Uangle);
    printf("Ukappa: %lf ", Ukappa);
    printf("Ubend: %lf ", stiff::Ubend);
    printf("Ucomp: %lf ", stiff::Ucomp);
    printf("Ushear: %lf ", stiff::Ushear);
    /* calc_stress(); */
    /* printf("Pscalar: %lf ", Pscalar); */

    if (mu != 0.0 || print_order_params == 1) {
      printf("U_ms: %lf lc: %lf", U_ms, lc_order_param);
      // cout << " doi: " << doi_lc << " gub_maxreal: " << gubbins_max_real
      cout << " gub_maxreal: " << gubbins_max_real
           << " gub_largest: " << gubbins_largest;
    }
    if ((polylen != 0.0 && sites_per_LC != 0) || print_order_params == 1) {
      cout << " poly_lc_order: " << poly_lc_order;
    }
    if (polylen != 0) cout << " aniso ratio " << anisotropic_ratio;
    printf("\n");
    fflush(stdout);
  }

}


void print_log(FILE *file) {
  if (myrank == 0) {
    fprintf(file, "%d %lf %lf %lf ", step, Ubond, Uangle, Ukappa);
    fprintf(file, "%lf %lf %lf ", stiff::Ubend, stiff::Ucomp, stiff::Ushear);

    if (mu != 0.0 || print_order_params == 1) {
      fprintf(file, "%lf %lf ", U_ms, lc_order_param);
      // fprintf( file, "%lf ", doi_lc) ;
      fprintf(file, "%lf %lf", gubbins_max_real, gubbins_largest);
    }

    if ((polylen != 0.0 && sites_per_LC != 0) || print_order_params == 1) {
      fprintf(file, " %lf", poly_lc_order);
    }
    fprintf(file, " %lf", anisotropic_ratio);
    fprintf(file, "\n");
    fflush(file);
  }
}

bool prepare_rg_variables(bool rg_done) {
  if (not polylen == 0.0) {
    if (not rg_done) {
      send_all_x_to_root();
      calculate_polymer_segments();
      rg_order();
      calc_scalar_parameters();
      R_parallel_sq = calculate_R_parallel();
      // Rg2 = calculate_Rg2();
      R_perpendicular_sq = (Rg2 - R_parallel_sq) / 2;
      anisotropic_ratio = sqrt(R_parallel_sq / R_perpendicular_sq);
    }
  }
  return true;
}
