import os
import sys
import numpy as np
from pylab import *
from math import *
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
#matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

aval = float(sys.argv[1])
bval = float(sys.argv[2])
cval = float(sys.argv[3])
tvals = np.arange(100); T = len(tvals);
epsilonvals = np.zeros(T);
for t in tvals:
	epsilonvals[t] = pow(aval*(1.0+1.0*t/bval),-cval);

plt.plot(tvals,epsilonvals/epsilonvals[0]) 
plt.show();