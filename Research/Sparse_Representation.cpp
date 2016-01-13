//
//  main.cpp
//  Research
//
//  Created by Christopher Chan on 1/6/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
// instructions: set working directory to research folder....
#include <iostream>
#include "Dependencies/cnpy.h"
#include"Voigt.hpp"
#include"Sparse_Representation.h"
#include<complex>
#include "Testing.hpp"
#include "/opt/local/include/armadillo" //use macports to install
#include <iomanip>

void Sparse_Representation::load_file () {
    cnpy::NpyArray arr = cnpy::npy_load(fname);
    int numGalaxy = arr.shape[0]-1; //for redshift
    int redshiftpos = arr.shape[0];
    int pdfsize = arr.shape[1];
    double* loaded_data = reinterpret_cast<double*>(arr.data);
    int curr = 0;
    for(int i = 0;i<numGalaxy;i++){
        vector<double> row;
        for(int i =0;i<pdfsize;i++){
            row.push_back(loaded_data[curr]);
            curr++;
        }
        pdfs.push_back(row);
    }
    //retrieve the last row of data for z... equivalent of z[-1]
    vector<double> z;
    for(int i =0;i<pdfsize;i++){
        z.push_back(loaded_data[curr]);
        curr++;
    }
    double dz = z[1]-z[0];
    arr.destruct();
    prepareDictionary(dz, numGalaxy, pdfsize,pdfs,z);
    
}
void Sparse_Representation::prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<vector<double> >& pdfs,vector<double>& z){
    //cout << dz;
    auto minelemz = std::min_element(std::begin(z),std::end(z));
    auto maxelemz = std::max_element(std::begin(z),std::end(z));
    tuple<double,double> mu = make_tuple(*minelemz,*maxelemz);
    int Nmu = 250;
    double max_sig = (*maxelemz-*minelemz)/12;
    double min_sig = dz/6;
    tuple<double,double> sig = make_tuple(min_sig,max_sig);
    int Nv = 3;
    int Nsig = 80;
    int Na = Nmu*Nsig*Nv;
    cout << "Nmu, Nsig, Nv = "<< "["<< Nmu<< ","<< Nsig<< ","<<Nv<< "]"<<endl;
    cout << "Total bases in dictionary: "<< Na<<endl;
    arma::Mat<double> D = create_voigt_dict(z, mu, Nmu, sig, Nsig, Nv);
    
    //Define some tolerance and or Max number of Basis to be used, when tolerance is reached rest of basis is 0.
    double toler = 1.e-10;
    int Nsparse = 20;
    int Ncoef = 32001;
    arma::vec AA = linspace(0,1,Ncoef);
    double Da = AA[1]-AA[0];
    
    //Python code does this in parallel. first doing this single threaded to ensure correctness'
    for(int ik = 0;ik<numGalaxy;ik++){
        arma::rowvec  pdf0 = D.row(ik);
        //cout << pdf0<<endl;
        break;
        
    }
    
}
void Sparse_Representation::sparse_basis(arma::Mat<double>& dictionary,arma::vec query_vec,int n_basis, int tolerance){
    
}
arma::Mat<double> Sparse_Representation::create_voigt_dict(vector<double>& zfine, tuple<double,double> mu, int Nmu, tuple<double,double> sigma, int Nsigma, int Nv,double cut){
    arma::vec zmid = linspace(get<0>(mu),get<1>(mu),Nmu);
    arma::vec sig = linspace(get<0>(sigma), get<1>(sigma), Nsigma);
    arma::vec gamma = linspace(0,.5,Nv);
    //cout<< setprecision(50)<< zmid.at(1)<<endl;
    int NA = Nmu*Nsigma*Nv;
    int Npdf = zfine.size();
    arma::Mat<double> A(Npdf,NA,arma::fill::zeros);
    int kk = 0;
    for(int i = 0; i<Nmu;i++){
        for(int j = 0;j<Nsigma;j++){
            for(int k = 0;k<Nv; k++){
                arma::vec pdft = Voigt::voigtprofile(zfine,zmid[i],sig[j],sig[j]*gamma[k]);
                //pdft = where(pdft>=cut,pdft,0.)
                for(int i = 0;i<pdft.size();i++){
                    if(pdft.at(i)<cut){
                        pdft(i) = 0;
                    }
                }
                double normalized = norm(pdft);
                arma::vec temp = pdft/normalized; //UNSURE ABOUT THIS PART OF CODE..CANT TEST CUZ DEBUGGER BAD. this could also be sped up dramatically....trying to ensure correctness first before speeding up code.
                for(int l = 0;l<temp.size();l++){
                    A(l,kk) = temp[l];
                }
                kk++;
                
            }
        }
    }
    return A;
}

double Sparse_Representation::norm(arma::vec& input){
    double accum = 0;
    for(int i = 0;i<input.size();i++){
        accum += input.at(i)*input.at(i);
    }
    return sqrt(accum);
}
//Helper functions
//http://stackoverflow.com/questions/11734322/matlab-type-arrays-in-c
arma::vec Sparse_Representation::linspace( double a,  double b, int n) {
    arma::vec array(n);
     double step = (b-a) / (n-1);
    int i = 0;
    while(a <= b) {
        array[i] = a;
        a += step;
        i++;// could recode to better handle rounding errors
    }
    return array;
}

int main () {
    Sparse_Representation init;
    init.load_file();
    Voigt temp;
    return 0;
}

