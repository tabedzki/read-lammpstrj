#include "globals.h"
#include "stiff-polymer-backbone.h"

void bond_setup(int, int, double, double);
void bond_setup(int, double, double);

void bonds_gauss() {
  int i, j, k, id1, id2, present_flag, m;
  double mdr2, mdr, dr[Dim], delr, mf;
  double ub_loc = 0.0, utmp;

  for (i = 0; i < ns_loc; i++) {
    id1 = my_inds[i];

    for (j = 0; j < n_bonds[id1]; j++) {
      id2 = bonded_to[id1][j];

      if (id2 < id1) continue;
      if (bond_coeff[id1][j] == 0) continue;
      mdr2 = pbc_mdr2(x[id1], x[id2], dr);
      mdr = sqrt(mdr2);
      delr = mdr - bond_eq[id1][j];
      ub_loc += delr * delr * bond_coeff[id1][j] / 2.0;
      mf = bond_coeff[id1][j] * delr / mdr;

      for (k = 0; k < Dim; k++) {
        f[id1][k] -= mf * dr[k];
        f[id2][k] += mf * dr[k];
      }
    }  // for ( j=0 ; j<n_bonds ;
  }    // for ( i=0 ; i<ns_loc

  // MPI_Allreduce sends the reduced value to all processors
  MPI_Allreduce(&ub_loc, &utmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Ubond = utmp;
}

/*
 * bonds_init sets up the bond list.
 *
 */

void bonds_init() {
  int i, j, k, m, n, ind = 0, id2, ind_old;

  for (i = 0; i < nstot; i++) n_bonds[i] = 0;

  for (j = 0; j < total_num_chains; j++) {
    if (do_gauss == 1) {
      for (m = 0; m < polylen - 1; m++) {
        ind = (j)*total_len + m;
        bond_setup(ind, poly_bond_req, poly_bond_k);
      }
    }

    for (m = 0; m < amount_of_side_chains; m++) {
      ind_old = (j)*total_len + (m)*freq_attached_group;

      if (spacer_len > 0) {
        ind = (j)*total_len + polylen + (m)*side_len;
        bond_setup(ind, ind_old, poly_spacer_bond_req, poly_spacer_bond_k);

        for (n = 0; n < spacer_len - 1; n++) {
          ind = (j)*total_len + polylen + (m)*side_len + n;
          bond_setup(ind, spacer_bond_req, spacer_bond_k);
        }
        ind_old = ind + 1;

      }  // end of spacer section

      ind = (j)*total_len + polylen + (m)*side_len + spacer_len;

      if (orient_head != 1) { ind = ind + (sites_per_LC - 1) / 2; }

      if (spacer_len > 0) {
        external_lc_bond_req = spacer_lc_bond_req;
        external_lc_bond_k = spacer_lc_bond_k;
      } else {
        external_lc_bond_req = poly_lc_bond_req;
        external_lc_bond_k = poly_lc_bond_k;
      }

      if (do_wlc == 1) { external_lc_bond_req = stiff::l0* stiff::gamma; }

      bond_setup(ind, ind_old, external_lc_bond_req, external_lc_bond_k);

      for (n = 0; n < sites_per_LC - 1; n++) {
        ind = (j)*total_len + polylen + (m)*side_len + spacer_len + n;

        bond_setup(ind, lc_bond_req, lc_bond_k);

      }  // end of lc_section
    }
  }  // end of total_num_chains

  for (i = 0; i < num_lc_only_c; i++) {
    for (m = 0; m < sites_per_LC - 1; m++) {
      ind = i * sites_per_LC + m + nsLCP;

      bond_setup(ind, lc_bond_req, lc_bond_k);
    }
  }
}

void bond_setup(int ind1, int ind2, double loc_bond_eq, double loc_bond_coeff) {
  bonded_to[ind1][n_bonds[ind1]] = ind2;
  bond_eq[ind1][n_bonds[ind1]] = loc_bond_eq;
  bond_coeff[ind1][n_bonds[ind1]] = loc_bond_coeff;
  n_bonds[ind1]++;

  bonded_to[ind2][n_bonds[ind2]] = ind1;
  bond_eq[ind2][n_bonds[ind2]] = loc_bond_eq;
  bond_coeff[ind2][n_bonds[ind2]] = loc_bond_coeff;
  n_bonds[ind2]++;
}

void bond_setup(int ind, double loc_bond_eq, double loc_bond_coeff) {
  bond_setup(ind, ind + 1, loc_bond_eq, loc_bond_coeff);
}