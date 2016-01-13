//
//  Sparse_Representation.h
//  Research
//
//  Created by Poop on 1/13/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#ifndef Sparse_Representation_h
#define Sparse_Representation_h
using namespace std;

class Sparse_Representation {
public:
    Sparse_Representation();
    void load_file();
    int area() {return 0;}//return width*height;}
    void prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<vector<double> >& pdfs,vector<double>& z);
    arma::Mat<double> create_voigt_dict(vector<double> &zfine, tuple<double,double> mu, int Nmu, tuple<double,double> sigma, int Nsigma, int Nv,double cut = 1.e-5);
    void sparse_basis(arma::Mat<double>& dictionary,arma::vec query_vec,int n_basis, int tolerance = 0);
private:
    string fname = "Data/CFHTLens_sample.P.npy";
    vector<vector<double> >pdfs;
    arma::vec linspace( double a,  double b, int n);
    double norm(arma::vec& input);

};
Sparse_Representation::Sparse_Representation(){
    //load_file();
}

#endif /* Sparse_Representation_h */
