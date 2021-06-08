#include <Eigen/Dense>

#ifdef STIFF

namespace stiff{
  double eta, gamma, epsilon_perp, epsilon_parallel, epsilon_b, l0;
  double Ubend, Ucomp, Ushear;
  double mdr_rij, mdr_rij2, mdr_rijrijjj_inv, temp_double, mdr_ririjj_inv; 
  Matrix<double, Dim, Dim> dRperp_drij, dRperp_drijj, dRperp_drijjj;
  Matrix<double, Dim, Dim> Ucross, temp_tensor;
  Matrix<double, Dim, 1>   Rij, Rperp, temp, temp2, uij, uijj, tf;
  Map<Matrix<double, Dim, 1>> ri(NULL), rij(NULL), rijj(NULL), rijjj(NULL);
  Map<Matrix<double, Dim, 1>> uij_raw(NULL), uijj_raw(NULL);
  
  void shear(double, int, Matrix<double, Dim, Dim>,
                    Matrix<double, Dim, 1>, Matrix<double, Dim, 1>);
  void calculate();
  void alt_compression();
  void alt_coupling();
  void alt_energy();

  namespace energy{
    double bending(int);
    double compression(int);
    double shear(int);
  };

  namespace forces{
  void shear(int, bool);
  void compression(int, bool);
  void bending(int, bool, bool);
  }
}

#else

namespace stiff{
  extern double eta, gamma, epsilon_perp, epsilon_parallel, epsilon_b, l0;
  extern double Ubend, Ucomp, Ushear;
  extern double mdr_rij, mdr_rij2, mdr_rijrijjj_inv, temp_double;
  extern Matrix<double, Dim, Dim> dRperp_drij, dRperp_drijj, dRperp_drijjj;
  extern Matrix<double, Dim, Dim> Ucross, temp_tensor;
  extern Matrix<double, Dim, 1>  Rperp, temp, temp2, uij, uijj,tf;
  // extern Map<Matrix<double, Dim, 1>> Rij, rij, rijj, rijjj;
  // extern Matrix<double, Dim, 1>  Rperp, temp, temp2;
  // extern Map<Matrix<double, Dim, 1>> Rij, rij, rijj, rijjj;
  extern Map<Matrix<double, Dim, 1>> uij_raw, uijj_raw;
  
  void shear(double, int, Matrix<double, Dim, Dim>,
                    Matrix<double, Dim, 1>, Matrix<double, Dim, 1>);
  void calculate();
  void alt_compression();
  void alt_coupling();
  void alt_energy();

  namespace energy{
    double bending(int);
    double compression(int);
    double shear(int);
  };

  namespace forces{
  void shear(int, bool);
  void compression(int, bool);
  void bending(int, bool, bool);
  }
}

#endif