//
//  LDASGMCMCmodel.h
//  
//
//  Created by Ayan Acharya on 9/20/15.
//
//

#ifndef ____LDASGMCMCmodel__
#define ____LDASGMCMCmodel__

#include "mathutils.h"
#include "samplers.h"
#include "LDASGMCMCdata.h"

class model
{
    unsigned int D,V,K,S,CollectionITER,BurninITER;
    double azero,bzero,czero,dzero,c,eta,gammazero,M;
    mat thetadk,thetatestdk,thetadkss,thetatestdkss,phiwk,phiwkss;
    rowvec thetakss,thetadss,cd,rk,rkss,xkss,Mk,ellkss,logpkss;
    mat xwkss,xdk;
    
public:
    
    model(unsigned int,unsigned int,double);
    colvec ProjSimplex(colvec); 
    void updatelocal(gsl_rng*,unsigned int,unsigned int,double,double,data);
    void updateglobal(gsl_rng*,double,double);
    void printresults(string,unsigned int);
    void seedmodel(string);
};

#endif /* defined(____LDASGMCMCmodel__) */