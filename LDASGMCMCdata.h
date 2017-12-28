//
//  rpfmdata.h
//  
//
//  Created by Ayan Acharya on 9/20/15.
//
//

#ifndef ____LDASGMCMCdata__
#define ____LDASGMCMCdata__

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <array>
#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

class data
{
	public:
		unsigned int D,Dtest,V,S;
		sp_mat Xdw,Xtestdw; 
		data(string);
};

#endif 
