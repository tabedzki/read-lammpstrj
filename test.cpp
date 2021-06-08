#include <Eigen/Dense>
// #include <eigen3/Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;

int main(void){

    int i, j, id, k,  n ;
    double  total=0;
    double l_total=0;
    EigenSolver<MatrixXf> ces ;
    ArrayXf Eigen_vals(3);
    ArrayXXf director;

    Eigen_vals(2)= 4.5;
    cout << Eigen_vals << endl;
    cout << Eigen_vals.cast<int>() << endl;
    cout << Eigen_vals.cast<int>().cast<double>() << endl;

    ArrayXd yo(3);
    yo[0] = 1;
    yo[1] = 2;
    yo[2] = 3;
    // yo << 1,4 ;
    yo *=2;
    cout << yo << endl;
    int ryan[3];
    ryan[0] = 5; 
    ryan[2] = 3; 
    ryan[1] = 9; 
    VectorXf yp(3);
    cout << *ryan << endl;
    for (int i =0; i < 3; i++) yp[i] = ryan[i]; 
    // yo(0) = ryan[0];
    yp /= 4;
    cout << yp  << endl;
    yp =  yp.array() * yp.array() ;
    // yp = yp.array() % 4;
    cout << yp  << endl;
    cout<< endl << endl;

    MatrixXf m(3,3);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
        m(i,j) =i; 
    }
    }
    cout << m << endl;
    cout <<m.sum() <<endl;
    
    cout << "Begin new region" << endl;
    ces.compute(m); 
    cout << "computed" << endl;
    Eigen_vals= m.eigenvalues().array().real().cast<float>();
    cout << "Eigenvals: " << endl;
    cout << Eigen_vals  << endl;
    double gubbins_max_real = Eigen_vals.maxCoeff(); 
    cout << "gubbins_max_real" << gubbins_max_real << endl;
    Eigen_vals.abs().maxCoeff(&n); 
    double gubbins_largest= Eigen_vals(n);
    cout << "bugger" << endl;
    director = ces.eigenvectors().array().real().cast<float>();
    cout << "n" << n << endl;
    cout << "gubbins_largest" << gubbins_largest << endl;
    cout << "pure eigenvectors:" << ces.eigenvectors() << endl;
    cout <<"director: " <<  director << endl;

    Eigen_vals = ces.eigenvalues().array().real().cast<float>();
    Eigen_vals.abs().maxCoeff(&n); 
    gubbins_largest= Eigen_vals(n);
    cout << "n" << n << endl;
    cout << "gubbins_largest" << gubbins_largest << endl;
    director = ces.eigenvectors().col(n).real();
    cout <<"director: " <<  director << endl;

    double **test = (double **) calloc(2,sizeof(double*));
    test[0] = (double *) calloc(6, sizeof(double));
    test[1] = (double *) calloc(6, sizeof(double));
    test[0][0] = 4;
    test[0][1]= 44;;
    test[0][2] = 55;
    test[0][3] = 66;
    test[0][4]= 11;
    test[0][5] = 99;

    test[1][0] = 1;
    test[1][1]= 3;
    test[1][2] = 9;
    test[1][3] = 7;
    test[1][4]= 11;
    test[1][5] = 8;
    /* cout << *test << endl; */
    cout << test[0][0] << endl;
    cout << test[0][1] << endl;
    cout << test[0][2] << endl;
    cout << test[0][3] << endl;
    cout << test[0][4] << endl;
    cout << test[0][5] << endl;
    cout << "hi"<< endl;

    /* yp( test,1,3); */
    Map<Matrix<double,1,3> > hehe  (test[0], 3);
    cout << hehe << endl;

    new (&hehe)  Map<Matrix<double,1,3> >  (test[1]+1, 3 );
    /* new (&hehe)  Map<Array<double,1,3> >  (test+1, 3 ); */
    cout << hehe << endl;
    cout << hehe.normalized() << endl;
    cout << hehe<< endl;
    cout << test[1][1] << endl;
    cout << test[1][2] << endl;
    cout << test[1][3] << endl;
    hehe.normalize();
    cout << hehe<< endl;
    cout << test[1][1] << endl;
    cout << test[1][2] << endl;
    cout << test[1][3] << endl;
return 0;
}
