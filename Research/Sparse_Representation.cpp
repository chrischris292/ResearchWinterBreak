//
//  main.cpp
//  Research
//
//  Created by Christopher Chan on 1/6/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
// instructions: set working directory to research folder....
#include"Sparse_Representation.h"

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
    Eigen::MatrixXd D = prepareDictionary(dz, numGalaxy, pdfsize,pdfs,z);
    
    //Define some tolerance and or Max number of Basis to be used, when tolerance is reached rest of basis is 0.
    double toler = 1.e-10;
    int Nsparse = 20;
    int Ncoef = 32001;
    Eigen::VectorXd AA;
    AA.setLinSpaced(Ncoef,0,1);
    double Da = AA[1]-AA[0];
    
    //Python code does this in parallel. first doing this single threaded to ensure correctness'
    for(int ik = 0;ik<numGalaxy;ik++){
        Eigen::VectorXd  pdf0(pdfs[ik].data(),pdfs[ik].size());
        int np = Nsparse;
        sparse_basis(D,pdf0,np);
        //cout << pdf0<<endl;
        break;
    }
    

    
}
/*
Eigen::MatrixXd Sparse_Representation::prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<vector<double> >& pdfs,vector<double>& z){
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
    Eigen::MatrixXd D = create_voigt_dict(z, mu, Nmu, sig, Nsig, Nv);
    return D;
    
}

void Sparse_Representation::sparse_basis(Eigen::MatrixXd& dictionary,Eigen::VectorXd query_vec,int n_basis, int tolerance){
    Eigen::VectorXd a_n(dictionary.rows());   //arma::zeros(dictionary.n_cols);
    Eigen::MatrixXd temp = dictionary.transpose();
    auto alpha = temp*query_vec.transpose();
    Eigen::MatrixXd alphaRes = alpha.eval();
    Eigen::VectorXd res = query_vec;
    Eigen::VectorXd idxs;
    idxs.setLinSpaced(dictionary.cols(),0,dictionary.cols());
    //cout << "idxs: "<< idxs<<endl;
}
Eigen::MatrixXd Sparse_Representation::create_voigt_dict(vector<double>& zfine, tuple<double,double> mu, int Nmu, tuple<double,double> sigma, int Nsigma, int Nv,double cut){
    Eigen::VectorXd zmid;
    zmid.setLinSpaced(Nmu,get<0>(mu),get<1>(mu));
    Eigen::VectorXd sig;
    sig.setLinSpaced(Nsigma,get<0>(sigma), get<1>(sigma));
    Eigen::VectorXd gamma;
    gamma.setLinSpaced(Nv,0,.5);
    //cout<< setprecision(50)<< zmid.at(1)<<endl;
    int NA = Nmu*Nsigma*Nv;
    int Npdf = zfine.size();
    Eigen::MatrixXd A(Npdf,NA);
    int kk = 0;
    for(int i = 0; i<Nmu;i++){
        for(int j = 0;j<Nsigma;j++){
            for(int k = 0;k<Nv; k++){
                Eigen::VectorXd pdft = Voigt::voigtprofile(zfine,zmid[i],sig[j],sig[j]*gamma[k]);
                //pdft = where(pdft>=cut,pdft,0.)
                for(int i = 0;i<pdft.size();i++){
                    if(pdft(i)<cut){
                        pdft(i) = 0;
                    }
                }
                double normalized = norm(pdft);
                Eigen::VectorXd temp = pdft/normalized; //UNSURE ABOUT THIS PART OF CODE..CANT TEST CUZ DEBUGGER BAD. this could also be sped up dramatically....trying to ensure correctness first before speeding up code.
                for(int l = 0;l<temp.size();l++){
                    A(l,kk) = temp[l];
                }
                kk++;
            }
        }
    }
    
    return A;
}

double Sparse_Representation::norm(Eigen::VectorXd& input){
    double accum = 0;
    for(int i = 0;i<input.size();i++){
        accum += input(i)*input(i);
    }
    return sqrt(accum);
}
//Helper functions
//http://stackoverflow.com/questions/11734322/matlab-type-arrays-in-c
Eigen::VectorXd Sparse_Representation::linspace( double a,  double b, int n) {
    Eigen::VectorXd array(n);
     double step = (b-a) / (n-1);
    int i = 0;
    while(a <= b) {
        array(i) = a;
        a += step;
        i++;// could recode to better handle rounding errors
    }
    return array;
}
*/
int main () {
    //Sparse_Representation init;
    ////init.load_file();
    //Voigt temp;
    return 0;
}

