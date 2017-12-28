// @ Ayan Acharya
// Date: 12/25/2017
// Code for LDA-SGMCMC


#include "LDASGMCMCdata.h"
#include "LDASGMCMCmodel.h"

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

	cout<<"LDA-SGMCMC main function called successfully.."<<endl;

	model LDASGMCMC(K,V,eta); batchnum = 0;
	for (directory_iterator itr(srcDirname); itr!=directory_iterator(); ++itr)
	{
	    trFilename = itr->path().string(); cout <<  trFilename << endl; cout<< (batchnum+1) <<endl;
		// data loading
		data Data(trFilename);
	    epsilont = pow(aval*(1.0 + (batchnum+1.0)/bval), -cval); // the learning rate at the t-th iteration
	    rhot     = 1.0;//1.0*Dapproxtotal/Data.D; 	             // scale by which the summary statistics from the t-th minitach needs to be upgraded
	    cout<<epsilont<<endl;
	    if (epsilont>=1)
	    {
	    	cout<<"learning rate is no smaller than 1, resetting it to 0.99"<<endl;
	    	epsilont = 0.99;
	    }
		// update local variables
		LDASGMCMC.updatelocal(rng, BurninITER, CollectionITER, epsilont, rhot, Data);
		// update global variables   
		LDASGMCMC.updateglobal(rng, epsilont, rhot);
		// print results for LDA-SGMCMC
		LDASGMCMC.printresults(opDirname, batchnum);
		batchnum = batchnum + 1;
	}
	//gsl_rng_free (rng);
	cout<<"program terminates .."<< endl;
	
	return 0;
}

