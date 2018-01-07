//
//  mathutils.h
//  
//
//  Created by Ayan Acharya on 9/16/15.
//
//

#ifndef ____mathutils__
#define ____mathutils__

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <armadillo>
#include <iterator>
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#define LOWLIMIT 1e-30
#define UPPERLIMIT 20

using namespace std;
using namespace arma;
using namespace boost;
using namespace boost::filesystem;

double logguard(double m);
double minguard(double m);

#endif /* defined(____mathutils__) */
