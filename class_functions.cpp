#include "globals.h"
#include "classes.hpp"

using namespace Eigen;
using namespace std; 

block::block(int type, int start_indx)
{
  type_molecule = type;
  start_ind = start_indx;
  length = types_of_lengths[type];
  end_ind = this->start_ind + this->length;
  internal_bond_length = 1;
  // internal_bond_strength =  ;
  // internal_bond_length = ;
  // angle_coeff = ;
}

block::block(void){
  die("Please initialize me with a type number.");
}

block::~block()
{
}

void block::main_chain_init(){
  int i;
  for (i = 0; i < Dim; i++){
    x[this->start_ind][i] = ran2() * L[i];
  }
};

void block::build_block(int N) 
{
  int i, ind, j;

  // Generate random_orientation vector
  double  mdu = 0.0 ; 
  VectorXd rand_u(Dim); 
  for ( j=0 ; j<Dim ; j++ ) {
    rand_u[j] = gasdev2() ;
  }
  rand_u.normalize(); 

  for(ind = 1 + this->start_ind; ind < this ->end_ind ; ind++){
    for( j = 0; j < Dim; j++){
      x[ind][j] = x[ind-1][j] + internal_bond_length * rand_u[j];
      // x[ind][j] = x[ind-1][j] + int_bond_length * rand_u[j]/mdu;        
      if ( x[ind][j] > L[j] ){
          x[ind][j] -= L[j] ;
      }
      else if ( x[ind][j] < 0.0 ){
          x[ind][j] += L[j] ;
      }
    } // for (j=0; j< Dim ...)
  } // for (ind = 1 ...)
  
};

double block::angle_calc(int ind1){
  int ind2 = ind1+1, ind3 = ind2+1;
  angle_calc(ind1, ind2, ind3) ;
}

double block::angle_calc(int ind1, int ind2, int ind3){
  VectorXd dr1(Dim);
  VectorXd dr2(Dim);
  vector <int> values = {ind1, ind2, ind3}; 
  double weight = calc_weight(values);
  if (weight ==0) return 0;
  int i;
  double sum;
  double dot_, cos_angle;
  for (i =0 ; i < Dim ; i++){
    dr1(i) = x[ind1][i] - x[ind2][i]; 
    dr2(i) = x[ind3][i] - x[ind2][i]; 
  }
  dot_ = dr1.dot(dr2);
  cos_angle = (dr1.dot(dr2)/(dr1.norm()*dr2.norm())); 
  sum =  this->angle_coeff * (1.0 + cos_angle);

  double fi, fk;
  for (i =0 ; i < Dim ; i++){
    fi  = this-> angle_coeff * 
      (cos_angle * dr1(i)/dr1.squaredNorm() - dr2(i) / dr2.norm() * dr1.norm()  );
    fk  = this-> angle_coeff * 
      (cos_angle * dr2(i)/dr2.squaredNorm() - dr1(i) / dr2.norm() * dr1.norm()  );
    
    f[ind1][i] += fi ;
    f[ind2][i] += -( fi + fk ) ;
    f[ind3][i] += fk ;
  }
  return sum;
}


// void block::angle_init(void){
// };



void block::attached_init(int ind_prev_block){
  
  int j; 
  // Generate random_orientation vector
  double  mdu = 0.0 ; // rand_u[Dim],
  VectorXd rand_u(Dim); 
  for ( j=0 ; j<Dim ; j++ ) {
    rand_u[j] = gasdev2() ;
  }
  rand_u.normalize(); 
  bonded_list details; 
  block *attached_group; 
  attached_group->start_ind = (how_many_allocated ++); 
  
  attached_group->end_ind = attached_group->start_ind + attached_group-> length; 

};

/**
 * @brief Calculates the value of the external bonds, calls for the program to
 * go through the internal bonds and then moves on to the next external bond
 * 
 * @return double  sum
 */
double block::ext_bonds(){ 
  int i; 
  double sum = 0.0;

  for ( i =0; i < this-> list_of_bonded.size(); i++){
    
    sum += bond_calculation(list_of_bonded[i]->bonded_from, 
    list_of_bonded[i]-> bonded_to, 
    list_of_bonded[i]->connection_type);

    sum += list_of_bonded[i]-> bonded_block -> int_bonds();
  }
  return sum;
};

/**
 * @brief Calculates the internal bonds of the system and then moves onto the
 * external bonds. Extenral bonds returns 0.0 if there is nothing attached. 
 * 
 * @return double  sum
 */
double block::int_bonds()
{
  double sum = 0.0;
  for (int i = start_ind; i < end_ind -1 ; i++){
  sum += bond_calculation(i , i + 1, 0) ;
  }
  sum += this-> ext_bonds();
  return sum;
};

VectorXd block::periodic_boundary(VectorXd dr){
  // Calculates the difference vector and pushes it to stay in the box
  for (int i =0 ; i < Dim ; i++){
    if ( dr[i] > Lh[i] )
      dr[i] -= L[i] ;

    else if ( dr[i] <= -Lh[i] )
      dr[i] += L[i] ;
    }
  return dr;
}

double block::bond_calculation(int ind1, int ind2, int int_ext_flag){
  vector <int> values {ind1, ind2};
  // double weight = calc_weight_bond(ind1, ind2);
  double weight = calc_weight(values);
  
  int i;
  VectorXd dr(Dim); 
  double ub_loc = 0.0;

  // Calculates the difference vector and pushes it to stay in the box
  for (i =0 ; i < Dim ; i++){
    dr(i) = x[ind1][i] - x[ind2][i]; 
    //   if ( dr[i] > Lh[i] )
    //   dr[i] -= L[i] ;

    // else if ( dr[i] <= -Lh[i] )
    //   dr[i] += L[i] ;
    }
  dr = periodic_boundary(dr);
    // dr = (dr.array() + Lengths_half.array())  %  Lengths

  double mag_dr = dr.norm();
  double delr = mag_dr - bond_dists[int_ext_flag];
  ub_loc += delr * delr - bond_coeffs[int_ext_flag]/2.0; 

  double mf = bond_coeffs[int_ext_flag] * delr / mag_dr;

  for ( i=0 ; i<Dim ; i++ ) {
    f[ind1][i] -= mf * dr[i] ;
    f[ind2][i] += mf * dr[i] ;
  }
  return ub_loc;
};



double block::bond_calculation(int ind1){
  int ind2;
  bond_calculation(ind1, ind2, 0 ) ; 
}

void block::init_internal_bond_length(){
  this->internal_bond_length = internal_bond_lengths[this->type_molecule];
}

void block::attach_polymer_block(int type, int attach_ind, int N ){
  
};


double block::calc_weight_bond(int ind1, int ind2){
  if (local_flag[ind1] ==1 and local_flag[ind2] ==1 )
    return 1.0;
  else if (local_flag[ind1] ==1 or local_flag[ind2] ==1 )
    return 0.5;
  else
    return 0.0;  
};


double block::calc_weight_angle(int ind1, int ind2){
  if (local_flag[ind1] ==1 and local_flag[ind2] ==1 )
    return 1.0;
  else if (local_flag[ind1] ==1 or local_flag[ind2] ==1 )
    return 0.5;
  else
    return 0.0;  
};

double block::calc_weight(vector<int> values){
  double weight = 0.0;
  for (int i = 0; i < values.size(); i ++)
    weight += local_flag[values[i]];
  return weight / values.size();
  // if (local_flag[ind1] ==1 and local_flag[ind2] ==1 )
  //   return 1.0;
  // else if (local_flag[ind1] ==1 or local_flag[ind2] ==1 )
  //   return 0.5;
  // else
  //   return 0.0;  
};


bool block::is_block_local(){
  int i;
  for ( i = start_ind; i < end_ind ; i++)
    if ( local_flag[i]) return true;
  
  for ( i =0; i < this-> list_of_bonded.size(); i++){
  if ( list_of_bonded[i]-> bonded_block -> is_block_local()) 
      return true;
  }
  return false;
}

bool block::do_we_send_north(){
  int i;
  for ( i = start_ind; i < end_ind  ; i++)
    if (x[i][Dim -1] < z_min + send_buff ) return true;
  
  for ( i =0; i < this-> list_of_bonded.size(); i++){
  if (list_of_bonded[i]-> bonded_block -> do_we_send_north()) 
      return true;
  }
  return false;
}

bool block::do_we_send_south(){
  int i;
  for ( i = start_ind; i < end_ind  ; i++)
    if (x[i][Dim -1] < z_min + send_buff ) return true;
  
  for ( i =0; i < this-> list_of_bonded.size(); i++){
  if (list_of_bonded[i]-> bonded_block -> do_we_send_south()) 
      return true;
  }
  return false;
}




void block::ms_forces(){
  if (this->type_molecule !=2) return;
  int center = (this -> start_ind + this->end_ind)/2;


}



// void block::old_school_build(){
//   int i, j; 
//   for (i=block->start_ind; i < block->end_ind; i++)
//   {}
// }


// void block::old_school_attach(int N ){


// }
