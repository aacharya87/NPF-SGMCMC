//
//  mathutils.cpp
//  
//
//  Created by Ayan Acharya on 9/16/15.
//
//

#include "mathutils.h"

double logguard(double m)
{
    // provides guard against log(0)
    return log(max(LOWLIMIT,m));
};

double minguard(double m)
{
    // provides guard against number lower than LOWLIMIT
    return max(LOWLIMIT,m);
};
