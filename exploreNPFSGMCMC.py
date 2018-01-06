#!/usr/bin/env python

########################################################################
## actual command:
## python -W ignore exploreNPFSGMCMC.py 20 200 200 64 48 10 1.0 0.7 0.7 0.7 0
## arguments: 
## 1st. K: maximum number of latent factors  
## 2nd: burnin: number of burn-in iterartions for Gibbs sampling
## 3rd: collection: number of iterations for collections of samples
## 4th: D: maximum number of documents per mini-batch (the actual number may be less)
## 5th: V: size of vocabulary
## 6th: N: number of mini-bacthes 
## 7th: eta: Dirichlet hyper-parameter for topics
## 8th: aval: variable used in the learning rate
## 9th: bval: variable used in the learning rate
## 10th: cval: variable used in the learning rate defined as (a+t/b)^(-c)
## 11th: seedind: flexibility to seed the learning from previous run (support will be added later)
#######################################################################

import os
import sys
import datetime
import numpy as np
from pylab import *
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
#matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import collections
from operator import itemgetter
##%matplotlib inline

def plotresults(opFileName,batchnum,M,thetadk,phiwk,rk):

	numL = 2; numW = 3;
	Assignment = thetadk.transpose(); Mest = np.dot(Assignment,phiwk);

	plt.figure((batchnum+1));
	sp = plt.subplot(numL,numW,1);
	plt.imshow(M)
	plt.title("original matrix")
	colorbar()

	sp = plt.subplot(numL,numW,2)	
	plt.imshow(Mest)
	plt.title("estimated matrix")
	colorbar()

	sp = plt.subplot(numL,numW,3)	
	plt.imshow(Assignment)
	plt.title("thetadk")
	colorbar()

	sp = plt.subplot(numL,numW,4)	
	plt.imshow(phiwk.transpose())
	plt.title("phiwk")
	colorbar()

	sp = plt.subplot(numL,numW,6)	
	plt.stem(rk)
	plt.title("rk")	

	plt.savefig(opFileName, dpi=200)

	#plt.show()	
	##wait = input('aa')

def createrandmat(D,V):
	
	M = np.random.random_integers(5, size=(D,V));
	M[np.nonzero(M<2)] = 0;
		
	return M;

def createsynthmat(D,V):
	
	mm = np.ones((D,V))
	mm = np.rint(mm*1)
	tmp1  = mm[0:D/4,0:V/4]
	tmp11 = np.zeros((D/4,3*V/4))
	tmp2  = mm[0:D/2,0:V/2]
	tmp21 = np.zeros((D/2,V/4))
	tmp3  = mm[0:D/4,0:V/4]
	tmp31 = np.zeros((D/4,3*V/4))
	M1 = np.concatenate((tmp1,tmp11),axis=1)
	M2 = np.concatenate((tmp21,tmp2,tmp21),axis=1)
	M3 = np.concatenate((tmp31,tmp3),axis=1)
	M  = np.concatenate((M1,M2,M3),axis=0)
	M  = 5.0*M;
		
	return M;

def writetoFile(fileName,M):
	
	f=open(fileName,'w'); idx = np.squeeze(np.nonzero(np.ravel(M, order='C'))); [D,V] = M.shape;
	writestring = str(D)+'\t'+str(V)+'\t'+str(len(idx))+'\n'; f.write(writestring)
	for indices in idx:
		i = int(np.floor(indices/V)); j = int((indices - i*V));
		writestring = str(i)+'\t'+str(j)+'\t'+str(5.0)+'\n'; f.write(writestring);
	f.close()	

K = int(sys.argv[1])
BurnIn = int(sys.argv[2])
Collection = int(sys.argv[3])
D = int(sys.argv[4])
V = int(sys.argv[5])
N = int(sys.argv[6])
eta = float(sys.argv[7])
aval = float(sys.argv[8])
bval = float(sys.argv[9])
cval = float(sys.argv[10])
seedind = int(sys.argv[11])

os.system('clear') 
os.system('rm -rf *.txt')
os.system('rm -rf NPFSGMCMC')

# for experiments with real data explicitly provide "trfile.txt", "predfile1.txt" and "predfile2.txt" 
execstring  = 'g++ -std=c++0x -o NPFSGMCMC NPFSGMCMC.cpp NPFSGMCMCmodel.cpp NPFSGMCMCdata.cpp samplers.cpp mathutils.cpp '
execstring += '-larmadillo -llapack -lblas -lboost_filesystem -lboost_system `gsl-config --cflags --libs`' 
print "compiling NPF-SGMCMC .."

## create few directories to store the input files and results
os.system(execstring); currdir = os.getcwd()+'/';
filesuffix = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
workDir    = currdir+'work_'+filesuffix; srcDir = workDir+'/srcDir'; opDir  = workDir+'/opDir';
os.mkdir(workDir); os.mkdir(srcDir); os.mkdir(opDir); Mlist = [];

## create the training data
Dapproxtotal = 0;
for n in np.arange(N):
	Dsize = D; ##np.random.randint(int(0.5*D), high=D);
	#M = createrandmat(Dsize,V);
	M = createsynthmat(Dsize,V); 
	Mlist.append(M);
	writetoFile(srcDir+'/trfile'+str(n+1),M);
	Dapproxtotal += Dsize;

execstring  = './NPFSGMCMC '+srcDir+' '+opDir+' '+str(K)+' '+str(V)+' '+str(BurnIn)+' '+str(Collection)+' '+str(eta)+' '+str(Dapproxtotal)+' ';
execstring += str(aval)+' '+str(bval)+' '+str(cval)+' '+str(seedind)+' '+str(seedind);
print "running NPF-SGMCMC .."
os.system(execstring)

print "reading results from NPF-SGMCMC .."

numbatches = int(0.25*len([name for name in os.listdir(opDir) if os.path.isfile(os.path.join(opDir, name))]))

for batchnum in np.arange(numbatches):	

	rk      = np.loadtxt(opDir+'/rk_'+str(batchnum)+'.txt');
	thetadk = np.loadtxt(opDir+'/thetadk_'+str(batchnum)+'.txt');
	phiwk   = np.loadtxt(opDir+'/phiwk_'+str(batchnum)+'.txt');

	#print thetadk.shape, phiwk.shape, rk.shape
	opFileName = opDir+'/'+str(batchnum+1)+'.png';
	plotresults(opFileName,batchnum,Mlist[batchnum],thetadk,phiwk,rk)

print "reading of results ends .."



