#include "globals.h"

void angle_setup(int, int, int, double);
void angle_setup(int, double);
double cross_vecs(double*, double*, double*);
double skew_times_vec(double*, double*, double*);
bool check_similarity(int, int, int);

void angles(void) {
  int i, j, k, i1, i2, i3;
  double rij[Dim], rkj[Dim], mdrij, mdrkj, cos_th, dot, u_loc, mdrij2, mdrkj2;
  double cross[Dim], temp[Dim];
  double sin_th, cross_mag2, skew_mag12, skew_mag32;
  double skewed32[Dim], skewed12[Dim];
  bool do_cos;

  Uangle = u_loc = 0.0;

  //  return ;

  for (i = 0; i < ns_loc; i++) {
    i1 = my_inds[i];

    for (j = 0; j < n_angles[i1]; j++) {
      if (i1 != angle_first[i1][j]) continue;
      if (angle_coeff[i1][j] == 0) continue;

      i2 = angle_mid[i1][j];
      i3 = angle_end[i1][j];

      // Generate the vectors rij, rkj
      mdrij2 = pbc_mdr2(x[i1], x[i2], rij);
      mdrkj2 = pbc_mdr2(x[i3], x[i2], rkj);

      
      do_cos = check_similarity(i1, i2, i3);
      if (do_cos){

        // Take dot product and determine the magnitude of rij, rkj
        dot = 0.0;

        for (k = 0; k < Dim; k++) dot += rij[k] * rkj[k];

        mdrij = sqrt(mdrij2);
        mdrkj = sqrt(mdrkj2);

        // cosine of angle between vectors
        cos_th = dot / mdrij / mdrkj;

        u_loc += angle_coeff[i1][j] * (1.0 + cos_th);

        for (k = 0; k < Dim; k++) {
          double fi = angle_coeff[i1][j] *
                      (cos_th * rij[k] / mdrij2 - rkj[k] / mdrij / mdrkj);
          double fk = angle_coeff[i1][j] *
                      (cos_th * rkj[k] / mdrkj2 - rij[k] / mdrij / mdrkj);
          f[i1][k] += fi;
          f[i2][k] += -(fi + fk);
          f[i3][k] += fk;
        }

      } else {
        cross_mag2 = cross_vecs(rij, rkj, cross);
        sin_th = sqrt(cross_mag2 / mdrij2 / mdrkj2);

        skew_mag12 = skew_times_vec(rij, cross, skewed12);
        skew_mag32 = skew_times_vec(rkj, cross, skewed32);

        for (k = 0; k < Dim; k++) {
          double fi = angle_coeff[i1][j] * sin_th *
                      ( skewed32[k] / cross_mag2 + rij[k] / mdrij2);
          double fk = angle_coeff[i1][j] * sin_th *
                      (-skewed12[k] / cross_mag2 + rkj[k] / mdrkj2);
          f[i1][k] += fi;
          f[i2][k] += -(fi + fk);
          f[i3][k] += fk;
        }
      }

    }  // for j < n_angles[i1]

  }  // for ( i < ns_loc

  MPI_Allreduce(&u_loc, &Uangle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void angles_alt(void) {
  int i, j, k, i1, i2, i3;
  double rij[Dim], rkj[Dim], mdrij, mdrkj, cos_th, dot, u_loc, mdrij2, mdrkj2;
  double cross[Dim], temp[Dim];
  double sin_th, cross_mag2, skew_mag12, skew_mag32;
  double skewed32[Dim], skewed12[Dim], angle_0;
  bool straight;

  Uangle = u_loc = 0.0;

  //  return ;

  for (i = 0; i < ns_loc; i++) {
    i1 = my_inds[i];

    for (j = 0; j < n_angles[i1]; j++) {
      if (i1 != angle_first[i1][j]) continue;
      if (angle_coeff[i1][j] == 0) continue;

      i2 = angle_mid[i1][j];
      i3 = angle_end[i1][j];

      // Generate the vectors rij, rkj
      mdrij2 = pbc_mdr2(x[i1], x[i2], rij);
      mdrkj2 = pbc_mdr2(x[i3], x[i2], rkj);

      
      straight = check_similarity(i1, i2, i3);
      straight ? angle_0 = -1: angle_0 = 0;

        // Take dot product and determine the magnitude of rij, rkj
        dot = 0.0;

        for (k = 0; k < Dim; k++) dot += rij[k] * rkj[k];

        mdrij = sqrt(mdrij2);
        mdrkj = sqrt(mdrkj2);

        // cosine of angle between vectors
        cos_th = dot / mdrij / mdrkj;

        double temp = (cos_th - angle_0);
        u_loc += angle_coeff[i1][j] * temp*temp;

        for (k = 0; k < Dim; k++) {
          double fi = angle_coeff[i1][j] * 2* temp *
                      (cos_th * rij[k] / mdrij2 - rkj[k] / mdrij / mdrkj);
          double fk = angle_coeff[i1][j] * 2* temp *
                      (cos_th * rkj[k] / mdrkj2 - rij[k] / mdrij / mdrkj);
          f[i1][k] += fi;
          f[i2][k] += -(fi + fk);
          f[i3][k] += fk;
        }


    }  // for j < n_angles[i1]

  }  // for ( i < ns_loc

  MPI_Allreduce(&u_loc, &Uangle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void angle_init(void) {
  int i, m, ind = 0, j, indf, indm, indb, ind_perp, ind_poly;
  if (myrank == 0) cout << "Entered angle_init" << endl;
  for (i = 0; i < nstot; i++) n_angles[i] = 0;

  // Angles for the LC chains attached to the polymers
  for (j = 0; j < total_num_chains; j++) {
    if (do_wlc == 1) {
      for (m = 0; m < polylen - 2; m++) {
        ind = (j)*total_len + m;
        angle_setup(ind, ind + 1, ind + 2, 1);
      }
    }
    for (m = 0; m < amount_of_side_chains; m++) {
      if (do_hard_angles_LC == 1) {
        ////////////////////////////////////////////////////
        // setting up the four 90 angles for the polymer  //
        ////////////////////////////////////////////////////

        //////////////////////////////////////
        // Polymer - LC/spacer angle setup  //
        //////////////////////////////////////

        ind_poly = j * total_len + m * freq_attached_group;
        // ind = j * total_len + polylen + m * side_len;
        if ((ind_poly % total_len) % polylen == polylen - 1) {
          indf = -1;
        } else {
          // indf = j * total_len + m * freq_attached_group + 1;
          indf = ind_poly + 1;
        }

        indm = ind_poly;
        if ((ind_poly % total_len) % polylen == 0) {
          indb = -1;
        } else {
          indb = ind_poly - 1;
        }

        if (spacer_len != 0){
          ind_perp = (j)*total_len + polylen + (m)*side_len;
          if (indb != -1) angle_setup(indb, indm, ind_perp, lc_angle_k);
          if (indf != -1) angle_setup(indf, indm, ind_perp, lc_angle_k);
        } else {
          if (orient_head == 0) { // Attached via the side
            if (sites_per_LC < 3) die("LC not long enough.");
            indm = (j)*total_len + polylen + (m)*side_len +
                   (sites_per_LC - 1) / 2 + spacer_len;

            if (indf != -1) angle_setup(indf, ind_poly, indm, lc_angle_k);
            if (indb != -1) angle_setup(indb, ind_poly, indm, lc_angle_k);
          } else { // Attached by the head
            if (sites_per_LC < 3) die("LC not long enough.");
            indm = (j)*total_len + polylen + (m)*side_len + spacer_len;

            if (indf != -1) angle_setup(indf, ind_poly, indm, lc_angle_k);
            if (indb != -1) angle_setup(indb, ind_poly, indm, lc_angle_k);
          }
        }


        ///////////////////////////////
        // LC - polymer angle setup  //
        ///////////////////////////////

        if (orient_head == 0) {
          indm = (j)*total_len + polylen + (m)*side_len +
                 (sites_per_LC - 1) / 2 + spacer_len;
          indf = indm + 1;
          indb = indm - 1;

          if (spacer_len == 0)
            // ind_perp = j * total_len + m * freq_attached_group;
            ind_perp = ind_poly;
          else
            ind_perp = j * total_len + m * side_len + polylen +
                       (spacer_len - 1);

          angle_setup(indf, indm, ind_perp, lc_angle_k);
          angle_setup(indb, indm, ind_perp, lc_angle_k);
        } else {
          indm = (j)*total_len + polylen + (m)*side_len + spacer_len;
          indf = indm + 1;

          if (spacer_len == 0)
            ind_perp = j * total_len + m * freq_attached_group;
          else
            ind_perp = j * total_len + m * freq_attached_group + polylen +
                       (spacer_len - 1);

          angle_setup(indf, indm, ind_perp, lc_angle_k);
        }
      }
      for (i = 0; i < sites_per_LC - 2; i++) {
        ind = (j)*total_len + polylen + (m)*side_len + i + spacer_len;
        angle_setup(ind, lc_angle_k);
      }
    }
  }
  // rod block angles //
  for (i = 0; i < num_lc_only_c; i++) {
    for (m = 0; m < sites_per_LC - 2; m++) {
      ind = i * sites_per_LC + m + nsLCP;
      angle_setup(ind, lc_angle_k);
    }
  }  // for ( i=0 ; i<nt[k]
}

void angle_setup(int ind1, int ind2, int ind3, double bond_coeff) {
  angle_first[ind1][n_angles[ind1]] = ind1;
  angle_mid[ind1][n_angles[ind1]] = ind2;
  angle_end[ind1][n_angles[ind1]] = ind3;
  angle_coeff[ind1][n_angles[ind1]] = bond_coeff;
  n_angles[ind1]++;

  angle_first[ind2][n_angles[ind2]] = ind1;
  angle_mid[ind2][n_angles[ind2]] = ind2;
  angle_end[ind2][n_angles[ind2]] = ind3;
  angle_coeff[ind2][n_angles[ind2]] = bond_coeff;
  n_angles[ind2]++;

  angle_first[ind3][n_angles[ind3]] = ind1;
  angle_mid[ind3][n_angles[ind3]] = ind2;
  angle_end[ind3][n_angles[ind3]] = ind3;
  angle_coeff[ind3][n_angles[ind3]] = bond_coeff;
  n_angles[ind3]++;
}

void angle_setup(int ind, double angle) {
  angle_setup(ind, ind + 1, ind + 2, angle);
}

double cross_vecs(double* v1, double* v2, double* out){
  double sum = 0;

  if (Dim == 3){
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  for (int i = 0; i < Dim; i++){
    sum += out[i] * out[i];
  }
  return sum;
}

double skew_times_vec(double* skewing_vec, double* v2, double* out){
  double sum = 0;
  out[0] = v2[1]*skewing_vec[2] - v2[2]*skewing_vec[1];
  out[1] = -v2[0]*skewing_vec[2] + v2[2]*skewing_vec[0];
  out[2] = skewing_vec[0]*v2[1] - skewing_vec[1]*v2[0];

  for (int i = 0; i < Dim; i++){
    sum += out[i] * out[i];
  }
  return sum;
}

bool check_similarity(int i1, int i2, int i3){
  int type1 = mol_type[i1];
  int type2 = mol_type[i2];
  int type3 = mol_type[i3];

  if (type1 ==1 ) type1 =0;
  if (type2 ==1 ) type2 =0;
  if (type3 ==1 ) type3 =0;

  if ((type1 == type2 && type2 != type3 ) || (type1 != type2 && type2 == type3 ))
  return false;
  else 
  return true; 
}