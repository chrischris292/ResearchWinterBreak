//
//  pdf_storage.cpp
//  Research
//
//  Created by Christopher Chan on 1/8/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#include <iostream>
#include<complex>

using namespace std;


class PDF_Storage {
    vector<vector<double> >pdfs;

public:
    PDF_Storage();
    void load_file();
    vector<vector<double>> create_voigt_dict(double zfine,double mu,int Nmu,double sigma,int, Nsigma, int Nv,double cut);

};
