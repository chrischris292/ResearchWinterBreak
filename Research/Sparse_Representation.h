//
//  Sparse_Representation.h
//  Research
//
//  Created by Poop on 1/13/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#ifndef Sparse_Representation_h
#define Sparse_Representation_h
#include <iostream>
#include "Dependencies/cnpy.h"
#include"Voigt.hpp"
#include<complex>
#include "Testing.hpp"
#include <iomanip>
#include "Dependencies/Eigen/Eigen"

using namespace std;

class Sparse_Representation {
public:
    Sparse_Representation();
    void load_file();
    int area() {return 0;}//return width*height;}
    Eigen::MatrixXd prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<double>& z);
    Eigen::MatrixXd create_voigt_dict(vector<double> &zfine, tuple<double,double> mu, int Nmu, tuple<double,double> sigma, int Nsigma, int Nv,double cut = 1.e-5);
    void sparse_basis(Eigen::MatrixXd& dictionary,Eigen::VectorXd query_vec,int n_basis, int tolerance = 0);
private:
    string fname = "Data/CFHTLens_sample.P.npy";
    vector<Eigen::VectorXd >pdfs;
    Eigen::VectorXd linspace( double a,  double b, int n);
    double norm(Eigen::VectorXd & input);
    int argMax(Eigen::VectorXd& input);


};
Sparse_Representation::Sparse_Representation(){
    //load_file();
}

#endif /* Sparse_Representation_h */
