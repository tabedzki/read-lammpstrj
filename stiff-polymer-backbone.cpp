#include <iomanip>
#include "globals.h"
#define STIFF
#include "stiff-polymer-backbone.h"

using namespace Eigen;

double pbc_mdr2( Matrix<double, Dim ,1 >& );
double pbc_mdr2( Matrix<double, 1, Dim >& );

extern Matrix<double,Dim,Dim> Eye;

void stiff::calculate(){
  int id_type, id, i;
  double Ubend_loc = 0, Ucomp_loc =0 , Ushear_loc = 0;
  bool include_4th;
  bool is_1st;

  // uijj is the main orientation vector we need. 
  // uij is the alt one.

  for ( i=0 ; i<ns_loc ; i++ ) {
    id = my_inds[i] ;
    
    // This does not account for homopolymers, will need to be fixed for that.
    if (mol_type[id] == 0 || mol_type[id] == 1 ){
      if ((id % total_len ) % polylen  == (polylen -1)
      // || (id % total_len ) % polylen  == 0 
      )
      continue;
    } else {
      continue;
    }

    new (&ri)    Map<Matrix<double, Dim ,1 >> (x[id+2], Dim);
    new (&rij)   Map<Matrix<double, Dim ,1 >> (x[id+1], Dim);
    new (&rijj)  Map<Matrix<double, Dim ,1 >> (x[id], Dim);
    if (not (id % total_len ) % polylen  == 0 ){
    new (&rijjj) Map<Matrix<double, Dim ,1 >> (x[id -1], Dim);
    }

    Rij  = rij - rijj; 
    pbc_mdr2(Rij);

    if ((id % total_len ) % polylen  == (polylen - 2)){
      uij = rij - rijj;
      include_4th = false;
    } else {
      uij = ri -  rijj; 
      include_4th = true;
    }

    if ((id % total_len ) % polylen  == 0 ){
      uijj= Rij;
      is_1st = true;
    } else {
      uijj= rij  - rijjj;
      is_1st = false;
    }

    pbc_mdr2(uij);
    pbc_mdr2(uijj);

    mdr_rijrijjj_inv = 1/uijj.norm();
    mdr_ririjj_inv = 1/uij.norm();

    uij.normalize(); 
    uijj.normalize(); 

    new (&uij_raw)  Map<Matrix<double, Dim ,1 >> (mono_u[id+1], Dim);
    new (&uijj_raw) Map<Matrix<double, Dim ,1 >> (mono_u[id], Dim);

    uij_raw  = uij;
    uijj_raw = uijj;


    Ucross = uijj * uijj.transpose();
    Rperp = Rij - Ucross * Rij ;

    temp_tensor = mdr_rijrijjj_inv *
      ( Rij.dot(uijj) * (Eye - Ucross)
      + ((Eye - Ucross) * Rij) * uijj.transpose());
    
    dRperp_drij = Eye - Ucross -  temp_tensor;

    dRperp_drijj = Ucross - Eye;
    dRperp_drijjj= temp_tensor;

    temp_tensor = uij * uij.transpose() ;

    // Force calculations
    forces::bending(id, include_4th, is_1st);
    forces::compression(id, is_1st);
    forces::shear(id, is_1st);

    // Energy calculations
    Ubend_loc  += energy::bending(id);
    Ucomp_loc  += energy::compression(id);
    Ushear_loc += energy::shear(id);
  }

  MPI_Allreduce(&Ubend_loc,  &Ubend,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Ucomp_loc,  &Ucomp,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Ushear_loc, &Ushear, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}




void stiff::forces::bending(int id, bool include_last, bool is_first) {
    temp = epsilon_b / stiff::l0 * (uij - uijj - stiff::eta * Rperp);

    // Takes into account whether or not affect the final particle position.
  if (include_last) {

    // Particle j + 1 from the derivation
    temp2 = (mdr_ririjj_inv * (Eye - temp_tensor)) * temp;
    for (int j = 0; j < Dim; j++) { f[id + 2][j] -= temp2(j); }

    // Particle j  from the derivation
    temp2 = -(mdr_rijrijjj_inv * (Eye - Ucross) + eta * dRperp_drij) * temp;
    for (int j = 0; j < Dim; j++) { f[id + 1][j] -= temp2(j); }
  } else {

    // Particle j  from the derivation
    temp2 = (mdr_ririjj_inv * (Eye - temp_tensor) -
               mdr_rijrijjj_inv * (Eye - Ucross) + eta * dRperp_drij) * temp;
    for (int j = 0; j < Dim; j++) { f[id + 1][j] -= temp2(j); }
  }

  // Particle j - 1 from the derivation
  temp2 = (mdr_ririjj_inv * (temp_tensor - Eye) - eta * dRperp_drijj) * temp;
  for (int j = 0; j < Dim; j++) {
    f[id][j] -= temp2(j); } 

    // Particle j - 2 from the derivation
    temp2= - (mdr_rijrijjj_inv * (Ucross - Eye) 
            + eta* dRperp_drijjj) * temp;
    if (is_first) {
      for (int j = 0; j < Dim; j++) {
        f[id][j] -= temp2(j);
      }
    } else {
      for (int j = 0; j < Dim; j++) {
        f[id - 1][j] -= temp2(j);
      }
    }
}




void stiff::forces::shear(int id, bool is_first){
    // Particle j from the derivation
    temp = epsilon_perp/stiff::l0 * dRperp_drij * Rperp;
    for (int j = 0; j < Dim; j++){ f[id+1][j] -= temp(j); } 

    // Particle j - 1 from the derivation
    temp = epsilon_perp/stiff::l0 * dRperp_drijj * Rperp;
    for (int j = 0; j < Dim; j++){ f[id][j] -= temp(j); } 

    // Particle j - 2 from the derivation
    temp = epsilon_perp/stiff::l0 * dRperp_drijjj * Rperp;
    if (is_first) {
      for (int j = 0; j < Dim; j++) {
        f[id][j] -= temp(j);
      }
    } else {
      for (int j = 0; j < Dim; j++) {
        f[id - 1][j] -= temp(j);
      }
    }
}

void stiff::forces::compression(int id, bool is_first){
  // Particle j from the derivation
  temp = epsilon_parallel / stiff::l0 *
         (Rij.dot(uijj) - stiff::l0 * stiff::gamma) *
         (mdr_rijrijjj_inv * (Eye - Ucross) * Rij + uijj);
  for (int j = 0; j < Dim; j++) { f[id + 1][j] -= temp(j); }

  // Particle j - 1 from the derivation
  temp = epsilon_parallel / stiff::l0 *
         (Rij.dot(uijj) - stiff::l0 * stiff::gamma) * (-uijj);
  for (int j = 0; j < Dim; j++) { f[id][j] -= temp(j); }

  // Particle j - 2 from the derivation
  temp = epsilon_parallel / stiff::l0 *
         (Rij.dot(uijj) - stiff::l0 * stiff::gamma) *
         (mdr_rijrijjj_inv * (Ucross - Eye) * Rij);
    if (is_first) {
      for (int j = 0; j < Dim; j++) {
        f[id][j] -= temp(j);
      }
    } else {
      for (int j = 0; j < Dim; j++) {
        f[id - 1][j] -= temp(j);
      }
    }
}

double stiff::energy::bending(int id) {
  temp = (uij - uijj - eta * Rperp);
  return epsilon_b * 0.5 / l0 * temp.dot(temp);
}

double stiff::energy::compression(int id) {
  temp_double = (Rij.dot(uijj) - l0 * stiff::gamma);
  return epsilon_parallel *0.5/l0 * temp_double * temp_double;
}

double stiff::energy::shear(int id) {
  return epsilon_perp * 0.5 /l0 * Rperp.dot(Rperp);
}