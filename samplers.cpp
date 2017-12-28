//
//  samplers.cpp
//  
//
//  Created by Ayan Acharya on 9/16/15.
//
//

#include "samplers.h"

double sampleCRT(gsl_rng *rng, const double m, const double gammazero)
{
    double sum = 0, bparam;
    for (int i=0; i<m; i++)
    {
        bparam = gammazero/(i + gammazero);
        if(gsl_rng_uniform (rng)<=bparam)
            sum = sum + 1;
    }
    return sum;
};

unsigned int TruncPoisson(gsl_rng *rng, double lambda)
{
    unsigned int m;
    double PMF = 1.0, prob;
    
    if(lambda>=1)
    {
        while(1)
        {
            m = gsl_ran_poisson (rng, lambda);
            if(m>0)
                break;
        }
    }
    else
    {
        m = 1;
        if(lambda<=1e-6)
            lambda = 1e-6;
        while(1)
        {
            prob = pow(lambda,m)*exp(-lambda)/(m*(1-exp(-lambda)));
            if (prob/PMF>gsl_rng_uniform (rng))
                break;        
            PMF = PMF-prob;
            m = m+1;
        }
    }
    return m;
};

/*unsigned int TruncBessel(gsl_rng *rng, double alpha)
{
	int j, nu, uval, modeval, countval;
	double prob, PMF;
	
	nu      = -1; countval = 0;
    modeval = floor((sqrt(pow(alpha,2) + pow(nu,2)) - nu)/2);
    if (modeval<1)
		modeval = 1;
    PMF     = cyl_bessel_i(nu, alpha)*exp(-alpha)*pow(alpha/2,-nu); 
    
    //std::cout<<modeval<<"\t"<<PMF<<"\t"<<cyl_bessel_i(nu, alpha)<<std::endl;

	for (j=0;PMF>0;j++)
	{
		uval = (modeval + j);
		prob = exp(uval*(2*log(alpha)-log(4))-lgamma(uval)-lgamma(uval+1) -alpha);
		if (prob/PMF >= gsl_rng_uniform (rng))
		{
			countval = uval;
			break;
		}
		else
			PMF -= prob;
		if(j>0 & j<modeval)
		{
			uval = (modeval - j);
			prob = exp(uval*(2*log(alpha)-log(4))-lgamma(uval)-lgamma(uval+1) -alpha);
			if (prob/PMF >= gsl_rng_uniform (rng))
			{
				countval = uval;
				break;
			}
			else
				PMF -= prob;
		}
	}
	return countval;
}*/

