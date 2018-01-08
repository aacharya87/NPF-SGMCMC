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
    azero = 1.0, bzero = 1.0, czero = 1.0, dzero = 1.0, gammazero = 1.0, c = 1.0; 
    K = Kval; V = Vval; eta = etaval;   
    
    // initialization of global variables
    xwkss = mat(V,K); xkss = rowvec(K); ellkss = rowvec(K); // the "ss" suffix just highlights that these terms are averaged over multiple samples, not derived from the final sample  
    phiwk = mat(V,K); phiwk.fill(1.0/V); M = 0.0; rk = rowvec(K); rk.fill(1.0*gammazero/K); rksum = sum(rk);
    // initialization of sufficient statistics
    phiwkss = mat(V,K,fill::zeros); rkss = rowvec(K,fill::zeros); Mk = rowvec(K,fill::zeros); // M and Mk get re-initialized afterwards

    cout<<"NPF-SGMCMC model initialization ends .."<<endl;
    return;
};

colvec model::ProjSimplex(colvec Phik, colvec oldphi) 
{
    colvec Phiknew = Phik + (1.0 - sum(Phik))*oldphi;
    double tmpsum = 0.0; 
    // put the values back to positive orthant and provide guard for extremely low values
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

void model::updatelocal(gsl_rng *rng, unsigned int citer, unsigned int biniter, unsigned int batchnum, double epsilont, double rhot, data Data)
{    
    cout<<"sampling of local variables for NPF-SGMCMC begins .."<<endl;
    CollectionITER = citer; BurninITER = biniter;  
    double param1,param2,val;
    double *tmpvecparam,*tmpvecvardouble,tmpsum;
    unsigned int *tmpvecvarint;
    unsigned int k,d,w;
    sp_mat::const_iterator start, end, it; sp_mat tmpsp;
 
    // size of the current mini-batch
    D = Data.D; 
    // initialization of local latent variables
    thetadk = mat(D,K); thetadk.fill(1.0/K); cd = rowvec(D); cd.fill(10.0); xdk = mat(D,K); 
    // initialization of sufficient statistics
	thetadkss = mat(D,K,fill::zeros); thetadss = rowvec(D);  
    // these have to be reset, as information only passes through Mk's from one mini-batch to the next 
    xkss.fill(0.0); xwkss.fill(0.0); ellkss.fill(0.0); logpdss = 0.0;

    cout<<"number of documents: "<<D<<endl;
    cout<<"number of words: "<<V<<endl;
    cout<<"number of topics: "<<K<<endl;    

    // Gibbs sampling iteration starts
    for (int i=0; i< (CollectionITER + BurninITER); i++)
    {
        if(i==0 || (i+1)%100 == 0)
            cout<< "Iteration: "<<(i+1)<<" of "<<(CollectionITER + BurninITER)<<", K: "<<K<<endl;        
        // reset few statistics first
        xdk.fill(0.0); 
        for (d=0; d<D; d++)
        {
            // sampling of latent counts; O(SK)
            tmpsp = Data.Xdw.row(d); start = tmpsp.begin(); end = tmpsp.end(); thetadss(d) = 0.0;
            for(it = start; it != end; ++it)
            {
                w = it.col(); val = (*it); tmpvecvarint = new unsigned int [K]; tmpvecparam = new double [K]; tmpsum = 0.0;
                for (k=0; k<K; k++)
                {
                    param1 = thetadk(d,k); param2 = phiwk(w,k); 
                    *(tmpvecparam+k) = param1*param2; tmpsum += *(tmpvecparam+k);
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
                        xkss(k)    += *(tmpvecvarint+k); 
                        xwkss(w,k) += *(tmpvecvarint+k); 
                    }
                }
                free(tmpvecvarint); free(tmpvecparam);
            }
            // sampling of thetadk O(DK)
            for (k=0; k<K; k++)
            {
                param1 = rk(k) + xdk(d,k); param2 = 1.0/(cd(d) + 1.0);
                thetadk(d,k) = minguard(gsl_ran_gamma(rng,param1,param2)); thetadss(d) += thetadk(d,k); 
                //update sufficient statistics for theta
                if (i>=BurninITER)
                {
                    thetadkss(d,k) += thetadk(d,k)/CollectionITER;
                    // sample CRT variables for update of rk's
                    ellkss(k)  += 1.0*sampleCRT(rng,xdk(d,k),rk(k)); 
                }
            }   
            // sampling of cd; O(D)
            param1 = (czero + rksum); param2 = 1.0/(dzero + thetadss(d));
            cd(d)  = minguard(gsl_ran_gamma(rng,param1,param2));   
            if (i>=BurninITER)
                logpdss += 1.0*logguard(1.0 + 1.0/cd(d));              
        }
    }

    // statistics to be used to update the global variables
    xkss   = rhot*xkss/CollectionITER; xwkss = rhot*xwkss/CollectionITER; 
    ellkss = rhot*ellkss/CollectionITER; logpdss = rhot*logpdss/CollectionITER;
    if (batchnum==0)
    {
        Mk = xkss; M = logpdss; 
    }
    else
    {
        Mk = (1.0-epsilont)*Mk + epsilont*xkss; 
        M  = (1.0-epsilont)*M  + epsilont*logpdss; 
    }

    cout<<"sampling of local variables for NPF-SGMCMC ends .."<<endl;
    return;
};

void model::updateglobal(gsl_rng *rng, double epsilont, double rhot)
{
    cout<<"sampling starts for global variables .."<< endl;
    double param1,param2;
    unsigned int k,w;
    colvec oldphi;
    rksum = 0.0;
    // sample global variables
	for (k=0; k<K; k++)
    {                    
        // sample phiwk
        oldphi = phiwk.col(k); param1 = 1.0*epsilont/Mk(k);
        for(w=0; w<V; w++)
        {
            param1     = epsilont*((xwkss(w,k) + eta) - phiwk(w,k)*(xkss(k) + eta*V))/Mk(k);
            param2     = pow(2*epsilont*phiwk(w,k)/Mk(k),0.5); // variance
            phiwk(w,k) = phiwk(w,k) + param1 + param2*gsl_ran_gaussian(rng, 1.0); 
        } 
        // project onto the Simplex 
        phiwk.col(k) = ProjSimplex(phiwk.col(k),oldphi);
        // sample rk
        param1 = epsilont*((ellkss(k) + 1.0*gammazero/K) - rk(k)*(logpdss + bzero))/M; 
        param2 = pow(2*epsilont*rk(k)/M,0.5); // variance
        rk(k)  = fabs(rk(k) + param1 + param2*gsl_ran_gaussian(rng, 1.0));
        rksum += rk(k);              
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
