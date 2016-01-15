//
//  Voigt.hpp
//  Research
//
//  Created by Christopher Chan on 1/8/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#ifndef Voigt_hpp
#define Voigt_hpp
#include<string>
#include<stdexcept>
#include<sstream>
#include<vector>
#include<cstdio>
#include<typeinfo>
#include<iostream>
#include<cassert>
#include<zlib.h>
#include<map>
#include <stdio.h>
#include<math.h>
#include<complex>
#include "Dependencies/Eigen/Eigen"

using namespace std;

class Voigt{
public:
    vector<vector<double> > create_voigt_dict(double zfine,double mu, int Nmu, double sigma, int Nsigma, int Nv, double cut);
    static Eigen::VectorXd voigtprofile(vector<double> x, double x_mean,double sigma, double gamma);
    
};
#endif /* Voigt_hpp */
