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
        Eigen::VectorXd row(pdfsize);
        for(int i =0;i<pdfsize;i++){
            row[i] = loaded_data[curr];
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
    Eigen::MatrixXd D = prepareDictionary(dz, numGalaxy, pdfsize,z);
    
    //Define some tolerance and or Max number of Basis to be used, when tolerance is reached rest of basis is 0.
    double toler = 1.e-10;
    int Nsparse = 20;
    int Ncoef = 32001;
    Eigen::VectorXd AA;
    AA.setLinSpaced(Ncoef,0,1);
    double Da = AA[1]-AA[0];
    
    //Python code does this in parallel. first doing this single threaded to ensure correctness'
    for(int ik = 0;ik<numGalaxy;ik++){
        Eigen::VectorXd  pdf0 = pdfs[ik];
        int np = Nsparse;
        sparse_basis(D,pdf0,np);
        break;
    }
    

    
}
/*
 ## preparing dictionary
 mu = [min(z), max(z)]
 Nmu = 250 #len(z)
 max_sig = (max(z) - min(z)) / 12.
 min_sig = dz / 6.
 sig = [min_sig, max_sig]
 Nv = 3
 Nsig = 80
 NA = Nmu * Nsig * Nv #final number of basis
 if rank == 0:
 print 'Nmu, Nsig, Nv = ', '[', Nmu, ',', Nsig, ',', Nv, ']'
 print 'Total bases in dictionary', NA
 D = create_voigt_dict(z, mu, Nmu, sig, Nsig, Nv)
 bigD = {}
 */
Eigen::MatrixXd Sparse_Representation::prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<double>& z){
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
/*
    Incomplete due to difficulty with array splicing/Cholesky solve.
 
 
 */
void Sparse_Representation::sparse_basis(Eigen::MatrixXd& dictionary,Eigen::VectorXd query_vec,int n_basis, int tolerance){
    Eigen::VectorXd a_n(dictionary.rows());   //arma::zeros(dictionary.n_cols);
    Eigen::MatrixXd temp = dictionary.transpose();
    Eigen::VectorXd alpha = temp*query_vec;
    Eigen::VectorXd res = query_vec;
    Eigen::MatrixXd alphaRes = alpha.eval();
    Eigen::VectorXd idxs;
    Eigen::VectorXd gamma;
    int n_active;
    idxs.setLinSpaced(dictionary.cols(),0,dictionary.cols()-1);
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n_basis,n_basis);
    L(0,0) = 1;
    for( n_active = 0;n_active<n_basis;n_active++){
        Eigen::VectorXd absVectorParam = dictionary.transpose()*res;
        int lam = argMax(absVectorParam); //can we improve this?
        cout << "Lam: "<<lam<<endl;

        if(n_active>0){
            //L[n_active, :n_active] = dot(dictionary[:, :n_active].T, dictionary[:, lam])
            cout << "A: "<<dictionary.leftCols(n_active)<<endl;
            cout << "B: "<<dictionary.col(lam)<< endl;
            Eigen::VectorXd result = dictionary.leftCols(n_active).transpose()*dictionary.col(lam);
            cout <<"result: "<<result<<endl;
            L.row(n_active).leftCols(n_active) = result;

            Eigen::VectorXd b= L.row(n_active).leftCols(n_active);
            cout << "TriangularA:"<<endl;
            cout<< L.topLeftCorner(n_active,n_active)<<endl;
            cout<<"TriangularB"<<endl;
            cout << b<<endl;
            cout<<"result:"<<endl;
            //L.row(n_active).leftCols(n_active) =L.topLeftCorner(n_active,n_active).triangularView<Eigen::Lower>().solve(L.row(n_active).leftCols(n_active).transpose()).transpose();
            L.row(n_active).leftCols(n_active) =L.topLeftCorner(n_active,n_active).triangularView<Eigen::Lower>().solve(b);
            //L.row(n_active).leftCols(n_active).norm()
            Eigen::VectorXd nextRow = L.row(n_active).leftCols(n_active);

            double v = normUnsquared(nextRow);
            //double v = pow(L.row(n_active).leftCols(n_active).norm(),2);
            L(n_active,n_active) = sqrt(1-v);
            //x = A.triangularView<Lower>().solve(b);
        }
        //dictionary[:, [n_active, lam]] = dictionary[:, [lam, n_active]]
        dictionary.col(lam).swap(dictionary.col(n_active));
        //alpha[[n_active, lam]] = alpha[[lam, n_active]]
        swapVectorVar(alpha,lam,n_active);
        //idxs[[n_active, lam]] = idxs[[lam, n_active]]
        swapVectorVar(idxs,lam,n_active);
        //gamma = sla.cho_solve((L[:n_active + 1, :n_active + 1], True), alpha[:n_active + 1], overwrite_b=False)x = A.llt() .solve(b));  // A sym. p.d.      #include <Eigen/Cholesky>
         gamma = L.topLeftCorner(n_active+1, n_active+1).llt().solve(alpha.head(n_active+1));
        //res = query_vec - dot(dictionary[:, :n_active + 1], gamma)
        res = query_vec - (dictionary.leftCols(n_active+1)*gamma);
    }

}
void Sparse_Representation::swapVectorVar(Eigen::VectorXd &input, int one, int two){
    double temp = input[one];
    input[one] = input[two];
    input[two] = temp;
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
                Eigen::VectorXd temp = pdft/normalized; //. this could also be sped up dramatically....trying to ensure correctness first before speeding up code.
                for(int l = 0;l<temp.size();l++){
                    A(l,kk) = temp[l];
                }
                kk++;
            }
        }
    }
    
    return A;
}
//Implement python argMax...assuming largest value is a single index...
int Sparse_Representation::argMax(const Eigen::VectorXd & input){
    int maxIndice = 0;
    for(int i = 0;i<input.size();i++){
        if(abs(input[maxIndice])<abs(input[i])){
            maxIndice = i;
        }
    }
    return maxIndice;
}

double Sparse_Representation::norm(Eigen::VectorXd& input){
    double accum = 0;
    for(int i = 0;i<input.size();i++){
        accum += abs(pow(input(i),2));
    }

    return sqrt(accum);
}
double Sparse_Representation::normUnsquared(Eigen::VectorXd& input){
    double accum = 0;
    for(int i = 0;i<input.size();i++){
        accum += abs(pow(input(i),2));
    }
    return accum;
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
int main () {
    Eigen::VectorXd c(7);
    c<< 0.556384  ,   0.358266   , -0.101661  ,   0.627309 ,  -0.175332 ,    0.357057 ,-3.35015e-17;
    Sparse_Representation init;
    //cout<<setprecision(50)<< sqrt(abs(1-init.normUnsquared(c)));

    init.load_file();
    return 0;
}

