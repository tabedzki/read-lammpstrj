#include <Eigen/Dense>
#include "globals.h"
// #include "class_functions.hpp"
#include <vector>
#include <map>

using namespace std ;
using namespace Eigen ;
  extern int how_many_allocated;
  extern VectorXf Lengths;
  extern VectorXf Lengths_half;
class monomer; class spacer; class lc_group; class block;


struct bonded_list{
  public: 
  int bonded_from; 
  int bonded_to; /// list of where the current block bonds to the next
  block* bonded_block;/// pointer to the bonded block
  int connection_type;
} ; 

class block
{
private:
  double bond_calculation(int, int);
  VectorXd periodic_boundary(VectorXd dr);
  double calc_weight_angle(int, int );
  double calc_weight_bond(int, int );
  double calc_weight(vector <int>);
public:
  int overall_poly_pos;
  int start_ind;
  int end_ind;
  int length;
  vector <bonded_list* > list_of_bonded ;
  int bonded_here ;
  int type_molecule; //0 = polymer_A; 1=spacer; 2=lc_group; 3 = polymer_B;
  int attaching_index = -1;
  bool any_monomers_here; 
  double angle_coeff; 
//   int monomers_on_this_processor;
//   double bond_coeffs[4], bond_dists[4];

  double internal_bond_strength;
//   double poly_spacer_bond_strength; 
//   double poly_lc_bond_strength; 
//   double spacer_lc_bond_strength; 

  double internal_bond_length; 
//   double poly_spacer_bond_length;this->internal_bond_length = internal_lengths[this->type_molecule]; 
//   double poly_lc_bond_length; 
//   double spacer_lc_bond_length; 
  
  
  bool is_block_local();
  bool do_we_send_north();
  bool do_we_send_south();
  void angle_init();
  void attach_block(int, int, int); 
  void build_block(int);
  double angle_calc(int); 
  double angle_calc(int,int, int); 
  void main_chain_init();
  void attached_init(int );
  void init_bond_properties();
  double ext_bonds();
  double int_bonds();
  double bond_calculation(int );
  void ms_forces();
  double bond_calculation(int, int, int );
  void attach_polymer_block(int, int, int );
  void init_internal_bond_length();
  block(int, int ); //Used for determining which type block is
  block(/* args */);
  ~block();
};






class polymer
{
private:
    /* data */
public:
    vector <monomer> poly_backbone; 
    polymer(int N);
    ~polymer();
};




class monomer
{
private:
    /* data */
public:
     
    // int pos[Dim]; 
    VectorXd pos; 
    VectorXd orientation; 
    monomer(/* args */);
    ~monomer();
};





class lc_group
{
private:
    /* data */
public:
    bool attached_by_head;
    VectorXd pos;
    // double orientation[Dim];
    VectorXd orientation; 
    lc_group(/* args */);
    ~lc_group();
    void init_orient();
};



class spacer
{
private:
    /* data */
public:
    VectorXd pos;
    VectorXd orientation; 
    spacer(/* args */);
    ~spacer();
};


