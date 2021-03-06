// @ Ayan Acharya
// Date: 12/25/2017
// Code for LDA-SGMCMC


#include "NPFSGMCMCdata.h"
#include "NPFSGMCMCmodel.h"

int main(int argc, char **argv)
{
	string trFilename,opDirname,srcDirname;
	int aRand, BurninITER, CollectionITER, K, V, D, seedval, batchnum, Dapproxtotal;
	double eta, epsilont, rhot, aval, bval, cval;
	const gsl_rng_type *Temp;
	
	// gsl set-up
	gsl_rng_env_setup(); Temp = gsl_rng_default; gsl_rng *rng = gsl_rng_alloc(Temp); srand (time(NULL));

	// provides different seeds for different runs
	aRand = rand() % 10 + 1; gsl_rng_set (rng, aRand);
	srcDirname = argv[1]; opDirname = argv[2]; K = atoi(argv[3]); V = atoi(argv[4]);
	BurninITER = atoi(argv[5]); CollectionITER = atoi(argv[6]); eta = atof(argv[7]); 
	Dapproxtotal = atoi(argv[8]); aval = atof(argv[9]); bval = atof(argv[10]); cval = atof(argv[11]); seedval = atoi(argv[12]);

	cout<<"NPF-SGMCMC main function called successfully.."<<endl;

	model NPFSGMCMC(K,V,eta); batchnum = 0;
	for (directory_iterator itr(srcDirname); itr!=directory_iterator(); ++itr)
	{
	    trFilename = itr->path().string(); 
		// data loading
		data Data(trFilename);
	    epsilont = pow(aval + (batchnum+1.0)/bval, -cval)/pow(aval + 1.0/bval, -cval);   // the learning rate at the t-th iteration
	    rhot     = 1.0*Dapproxtotal/Data.D; 	             // scale by which the summary statistics from the t-th minitach needs to be upgraded
	    cout<< "filename: " << trFilename << endl; 
	    cout<< "batch number: "<<(batchnum+1) <<endl;
	    cout<< "learning rate: " << epsilont << endl;
	    cout<< "rhot: " << rhot << endl;
	    if (epsilont>=1)
	    {
	    	cout<<"learning rate is no smaller than 1, resetting it to 0.99"<<endl;
	    	epsilont = 0.99;
	    }
		// update local variables
		NPFSGMCMC.updatelocal(rng, BurninITER, CollectionITER, batchnum, epsilont, rhot, Data);
		// update global variables   
		NPFSGMCMC.updateglobal(rng, epsilont, rhot);
		// print results for LDA-SGMCMC
		NPFSGMCMC.printresults(opDirname, batchnum);
		batchnum = batchnum + 1;
	}
	//gsl_rng_free (rng);
	cout<<"program terminates .."<< endl;
	
	return 0;
}

