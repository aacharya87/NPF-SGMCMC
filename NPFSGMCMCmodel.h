//
//  NPFSGMCMCmodel.h
//  
//
//  Created by Ayan Acharya on 9/20/15.
//
//

#ifndef ____NPFSGMCMCmodel__
#define ____NPFSGMCMCmodel__

#include "mathutils.h"
#include "samplers.h"
#include "NPFSGMCMCdata.h"

class model
{
    unsigned int D,V,K,S,CollectionITER,BurninITER;
    double azero,bzero,czero,dzero,c,eta,gammazero,M,rksum,logpdss;
    mat thetadk,thetatestdk,thetadkss,thetatestdkss,phiwk,phiwkss;
    rowvec thetakss,thetadss,cd,rk,rkss,xkss,Mk,ellkss;
    mat xwkss,xdk;
    
public:
    
    model(unsigned int,unsigned int,double);
    colvec ProjSimplex(colvec,colvec); 
    void updatelocal(gsl_rng*,unsigned int,unsigned int,unsigned int,double,double,data);
    void updateglobal(gsl_rng*,double,double);
    void printresults(string,unsigned int);
    void seedmodel(string);
};

#endif /* defined(____NPFSGMCMCmodel__) */
