//
//  bigDTuple.hpp
//  Research
//
//  Created by Christopher Chan on 2/18/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#ifndef bigDTuple_hpp
#define bigDTuple_hpp

#include <stdio.h>;
#include "Dependencies/Eigen/Eigen";
using namespace std;
class bigDTuple{
public:
    tuple<Eigen::VectorXd,Eigen::VectorXd> sparse;
    Eigen::VectorXd sparse_ind;
};
#endif /* bigDTuple_hpp */
