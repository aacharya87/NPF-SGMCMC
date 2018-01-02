//
//  NPFSGMCMCmodel.cpp
//  
//
//  Created by Ayan Acharya on 12/25/17.
//
//

#include "NPFSGMCMCmodel.h"

model::model(unsigned int Kval, unsigned int Vval, double etaval)
{    
    cout<<"intializing NPF-SGMCMC model .."<<endl;
    
    // model hyper-parameters initialized
    azero = 1.0, bzero = 1.0, czero = 1.0, dzero = 1.0, gammazero = 1.0*azero/bzero, c = 1.0; 
    K = Kval; V = Vval; eta = etaval;   
    
    // initialization of global variables
    xwkss = mat(V,K); xkss = rowvec(K); ellkss = rowvec(K); logpkss = rowvec(K); // the "ss" suffix just highlights that these are averaged over multiple samples, not derived from the final sample  
    phiwk = mat(V,K); phiwk.fill(1.0/V); M = 0.0; rk = rowvec(K); rk.fill(1.0*gammazero/K); 
    // initialization of sufficient statistics
    phiwkss = mat(V,K,fill::zeros); rkss = rowvec(K,fill::zeros); Mk = rowvec(K,fill::zeros); 

    cout<<"NPF-SGMCMC model initialization ends .."<<endl;
    return;
};

colvec model::ProjSimplex(colvec Phik, colvec oldphi) 
{
    colvec Phiknew = Phik + (1.0 - sum(Phik))*oldphi;
    double tmpsum = 0.0; 
    // put the values back to positive orthant
    for (int w=0; w<V; w++)
    {
        if (Phiknew(w)<LOWLIMIT)
            Phiknew(w) = LOWLIMIT;
        tmpsum += Phiknew(w);
    }
    // normalize
    Phiknew = Phiknew/tmpsum;
    return Phiknew;
}

void model::updatelocal(gsl_rng *rng, unsigned int citer, unsigned int biniter, double epsilont, double rhot, data Data)
{    
    cout<<"sampling of local variables for NPF-SGMCMC begins .."<<endl;
    CollectionITER = citer; BurninITER = biniter;  
    double param1,param2,val,rksum;
    double *tmpvecparam,*tmpvecvardouble,tmpsum;
    unsigned int *tmpvecvarint;
    unsigned int k,d,w;
    sp_mat::const_iterator start, end, it; sp_mat tmpsp;
 
    // size of the current mini-batch
    D = Data.D; 
    // initialization of local latent variables
    thetadk = mat(D,K); thetadk.fill(1.0/K); cd = rowvec(D,fill::ones); xdk = mat(D,K); 
    // initialization of sufficient statistics
	thetadkss = mat(D,K,fill::zeros); thetadss = rowvec(D); rksum = sum(rk);
    // these have to be reset, as information only passes through Mk's from one mini-batch to the next 
    xkss.fill(0); xwkss.fill(0); ellkss.fill(0.0); logpkss.fill(0.0);

    cout<<"number of documents: "<<D<<endl;
    cout<<"number of words: "<<V<<endl;
    cout<<"number of topics: "<<K<<endl;    

    // Gibbs sampling iteration starts
    for (int i=0; i< (CollectionITER + BurninITER); i++)
    {
        if(i==0 || (i+1)%100 == 0)
            cout<< "Iteration: "<<(i+1)<<" of "<<(CollectionITER + BurninITER)<<", K: "<<K<<endl;        
        // reset few statistics first
        thetadss.fill(0.0); xdk.fill(0); 
        // sampling of latent counts; O(SK)
        for (d=0; d<D; d++)
        {
            tmpsp = Data.Xdw.row(d); start = tmpsp.begin(); end = tmpsp.end();
            for(it = start; it != end; ++it)
            {
                w = it.col(); val = (*it); tmpvecvarint = new unsigned int [K]; tmpvecparam = new double [K]; tmpsum = 0.0;
                for (k=0; k<K; k++)
                {
                    param1 = thetadk(d,k); param2 = phiwk(w,k); *(tmpvecparam+k) = param1*param2; tmpsum += *(tmpvecparam+k);
                }
                //normalization
                for (k=0; k<K; k++)
                    *(tmpvecparam+k) = *(tmpvecparam+k)/tmpsum;
                gsl_ran_multinomial(rng, K, val, tmpvecparam, tmpvecvarint);
                // update sufficient statistics of latent rates
                for (k=0; k<K; k++)
                {
                    xdk(d,k) += *(tmpvecvarint+k); 
                    if (i>=BurninITER)
                    {
                        xkss(k) += *(tmpvecvarint+k); xwkss(w,k) += *(tmpvecvarint+k); 
                    }
                }
                free(tmpvecvarint); free(tmpvecparam);
            }
        }
        // sampling of thetadk O(DK)
        for (k=0; k<K; k++)
        {
            for (d=0; d<D; d++)
            {
                // sample thetadk
                param1       = rk(k) + xdk(d,k); param2 = 1.0/(cd(d) + 1.0);
                thetadk(d,k) = minguard(gsl_ran_gamma(rng,param1,param2));
                //update sufficient statistics for theta
                thetadss(d) += thetadk(d,k);
                if (i>=BurninITER)
                {
                    thetadkss(d,k) += thetadk(d,k)/CollectionITER;
                    // sample CRT variables for update of rk's
                    ellkss(k) += sampleCRT(rng,xdk(d,k),rk(k)); logpkss(k) += logguard(1.0 + 1.0/cd(d));   
                }
            }   
        }   
        // sample cd; O(D)
        for (d=0; d<D; d++)
        {
            param1 = (czero + rksum); param2 = 1.0/(dzero + thetadss(d));
            cd(d)  = minguard(gsl_ran_gamma(rng,param1,param2));
        }
    }
    // statistics to be used to update the global variables
    xkss = xkss/BurninITER; xwkss = xwkss/BurninITER; ellkss = ellkss/BurninITER; logpkss = logpkss/BurninITER;
    Mk   = (1.0-epsilont)*Mk + epsilont*rhot*xkss;
    M    = (1.0-epsilont)*M  + epsilont*rhot*sum(ellkss); 

    cout<<"sampling of local variables for NPF-SGMCMC ends .."<<endl;
    return;
};

void model::updateglobal(gsl_rng *rng, double epsilont, double rhot)
{
    cout<<"sampling starts for global variables .."<< endl;
    double param1,param2,param3,param4;
    unsigned int k,w;
    colvec oldphi;
    
    // sample global variables
	for (k=0; k<K; k++)
    {                    
        // sample phiwk
        for(w=0; w<V; w++)
        {
            param1 = 1.0*(epsilont/Mk(k)); param2 = param1*(rhot*xwkss(w,k) + eta);
            param3 = (1.0 - param1*(rhot*xkss(k) + eta*V)); param4 = pow(2*param1*phiwk(w,k),0.5); 
            oldphi = phiwk.col(k);
            phiwk(w,k) = param2 + param3*phiwk(w,k) + gsl_ran_gaussian(rng, param4); 
        } 
        // project onto the Simplex 
        phiwk.col(k) = ProjSimplex(phiwk.col(k),oldphi);
        // sample rk
        param1 = 1.0*(epsilont/M); param2 = param1*(rhot*ellkss(k) + 1.0*gammazero/K); 
        param3 = (1.0 - param1*(rhot*logpkss(k) + bzero)); param4 = pow(2*param1*rk(k),0.5); 
        rk(k)  = fabs(param2 + param3*rk(k) + gsl_ran_gaussian(rng, param4));                
    }
    phiwkss += phiwk; rkss += rk;
    cout<<"sampling terminates for global variables .."<< endl;;
};

void model::printresults(string opDirname, unsigned int batchnum)
{
    int k,d,w;
    
    cout<<"printing results .."<< endl; cout.precision(10);
    
    ofstream opfile1(opDirname+"/thetadk_"+to_string(batchnum)+".txt"); 
    ofstream opfile2(opDirname+"/phiwk_"+to_string(batchnum)+".txt");		
    ofstream opfile3(opDirname+"/rk_"+to_string(batchnum)+".txt");     

    for (k=0; k<K; k++)
    {
        for(d=0;d<D;d++)
            opfile1<<thetadkss(d,k)<<"\t";
        opfile1<<endl;
        for(w=0;w<V;w++)
            opfile2<<phiwkss(w,k)/(batchnum+1)<<"\t"; // the normalization is necessary
        opfile2<<endl;   
        opfile3<<rkss(k)/(batchnum+1)<<"\t";     // the normalization is necessary
    }	
	opfile1.close(); opfile2.close(); opfile3.close();
	
	ofstream opfile4(opDirname+"/prediction_"+to_string(batchnum)+".txt"); 	
	/*for(d=0;d<D;d++)
	{
		for(w=0;w<V;w++)
			opfile4<<thetaphiss(d,w)<<"\t";
		opfile4<<endl;	
    }*/
    opfile4.close(); 
    
    cout<<"printing of results done.."<< endl;	
};
