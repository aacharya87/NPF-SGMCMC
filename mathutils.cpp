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
    if(m<LOWLIMIT)
        return log(LOWLIMIT);
    else
        return log(m);
};

double minguard(double m)
{
    // provides guard against number lower than LOWLIMIT
    if(m<LOWLIMIT)
        return (LOWLIMIT);
    else
        return m;
};
