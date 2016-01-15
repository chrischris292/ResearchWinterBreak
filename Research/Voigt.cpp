//
//  Voigt.cpp
//  Research
//
//  Created by Christopher Chan on 1/8/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#include "Voigt.hpp"
#include "Dependencies/Faddeeva.hpp"
using namespace std;

/*
Eigen::VectorXd Voigt::voigtprofile(vector<double> x, double x_mean,double sigma, double gamma){
    //x=x-x_mean
    for(int i = 0;i<x.size();i++){
        x[i] = x[i]-x_mean;
    }
    vector<complex<double>> zTemp;
    Eigen::VectorXd returnVal(x.size());
    for(int i = 0;i<x.size();i++){
        complex<double> temp(x[i],gamma);
        temp.operator/=(sqrt(2.)*sigma);
        zTemp.push_back(temp);
        returnVal[i] = Faddeeva::w(zTemp[i]).real();
    }
    return returnVal;
}
*/