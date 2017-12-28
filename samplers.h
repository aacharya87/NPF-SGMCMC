//
//  samplers.h
//  
//
//  Created by Ayan Acharya on 9/16/15.
//
//

#ifndef ____samplers__
#define ____samplers__

#include <math.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <boost/math/special_functions/bessel.hpp>
//using namespace boost::math;


double sampleCRT(gsl_rng *rng, const double m, const double gammazero);
unsigned int TruncPoisson(gsl_rng *rng, double lambda);
//unsigned int TruncBessel(gsl_rng *rng, double alpha);

#endif /* defined(____samplers__) */

