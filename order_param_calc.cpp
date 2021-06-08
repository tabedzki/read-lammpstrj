#include "globals.h"
double P2(double);
double sp_plain, sp_star, sq;

// Make sure that calculate_polymer_segments(); is invoked before this
void calc_scalar_parameters() {
  if (myrank == 0) {
    int id, mod_id, plain_cnt = 0, star_cnt = 0;
    sp_plain = 0, sp_star = 0, sq = 0;
    Map<Matrix<double, 1, Dim>> local_dr_store(NULL);
    Matrix <double,1,Dim> average_dr_store;
    int freq_attached_group_loc;

    new (&local_dr_store) Map<Matrix<double, 1, Dim>>(dr_store[0], Dim); //dr_stored stored
    for (int j = 0; j < total_num_chains; j++) {
      for (int m = 0; m < (polylen - 2); m++) {
      id = (j)*total_len + m;
        mod_id = id % total_len;
        // Allows for usage of Eigen commands on top of the C arrays.

      // Calculates 6a from Wang and Wang
        if (freq_attached_group != 0) {
          if ((mod_id < Nda) && (mod_id % freq_attached_group) == 0) {
            average_dr_store = local_dr_store;
            new (&local_dr_store) Map<Matrix<double, 1, Dim>>(dr_store[id+1], Dim); // dr_stored
            average_dr_store = average_dr_store + local_dr_store;
            sp_plain += P2(average_dr_store.normalized().dot(director));
        plain_cnt++;
          } else
      // Calculates 6b from Wang and Wang
              if (mod_id < polylen - 1) {
        sp_star += P2(local_dr_store.normalized().dot(director));
        star_cnt++;
          } else {
            die("Not part of the backbone!", 44);
      } 
        } else if (mod_id < polylen - 1) {
          sp_star += P2(local_dr_store.normalized().dot(director));
          star_cnt++;
        } else {
        die("Not part of the backbone!", 44);
      }
    }
    sp_plain /= plain_cnt;
    sp_star  /= star_cnt ;
  }
}
  MPI_Barrier(MPI_COMM_WORLD);
}

// Make sure that calculate_polymer_segments(); is invoked before this
double calculate_R_parallel(void) {
  if (myrank == 0) {
    double sum1 = 0;
    double Rg_parallel_calc = 0;
    int id;
    double ddot = 0;
    Map<Matrix<double, 1, Dim>> local_dr_store(NULL);

    for (int j = 0; j < total_num_chains; j++) {
      for (int m = 0; m < (polylen ); m++) {
        for (int n = m; n < (polylen); n++) {
          sum1 = 0;
          for (int o = m ; o < n; o++) {
            id = (j)*total_len + o;
            new (&local_dr_store) Map<Matrix<double, 1, Dim>>(dr_store[id]); //dr_stored

            ddot = local_dr_store.dot(director);
            sum1 += ddot;
          }
          Rg_parallel_calc += sum1 * sum1;
        }
      }
    }

    Rg_parallel_calc = Rg_parallel_calc / polylen / polylen / total_num_chains;

    return Rg_parallel_calc;
  }
  return 0;
}