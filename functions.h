#include "fftw3-mpi.h"
#include <Eigen/Dense>
#include <complex>

using namespace std ;

double ran2(void ) ;
double gasdev2( void ) ;
int cell_stack( int* );
void cell_unstack( int , int* );
void field_gradient( double* , double* , int ) ;
void field_gradient_cdif( double* , double* , int ) ;
void convolve_fields( double*, double*, double* ) ;
void write_grid_data( const char* , double* ) ;
void write_Sfield_data( const char* , double*** ) ;
void write_kspace_data( const char* , complex<double>* ) ;

FILE *fopen_mkdir(char*, char*);
void rek_mkdir(char*);

int unstack_stack( int ) ;
void unstack_local( int, int* ) ;
int stack( int* ) ;
int stack_to_local( int* ) ;
int stack_local( int* ) ;
double integrate( double* );
void unstack( int , int* );
double get_k( int , double* ) ;
double get_k_alias( int , double* ) ;
void get_r( int , double* ) ;

void fftw_fwd( double* , complex<double>* );
void fftw_back( complex<double>* , double* );
int remove_dupes( int* , int );

int malloc2ddouble( double***, int , int  ) ;
int malloc2dint( int***, int , int  ) ;

void pbc_vdr( double*, double* , double* );
double pbc_mdr2( double*, double* , double* );
void die(const char *kill);
void die(const char *kill, int);
void die();

int mpi_Bcast_posits( int, int, double**, int, MPI_Comm ) ;
int allocate_mpi_tmp( int, int ) ;
void find_passengers( void ) ;
int is_partic_in_list( int, int, int* ) ;

void add_to_array( int, int*, double**, double** ) ;
void insert_in_array( int, int*, double**, double** ) ;
void extract_from_array( int, int*, double**, double** ) ;

double double_dot( double**, double** ) ;
int LC_has_orientation( int ) ;
double  diag_mat( double**, double* ) ;
void normalize_vector(double *);
void zero_vector(double *);
double P2(double);
