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

aval = 0.1; bval = aval; cval = aval;
tvals = np.arange(10); T = len(tvals);
epsilonvals = np.zeros(T);
for t in tvals:
	epsilonvals[t] = pow(aval*(1+t/bval),-cval);

plt.plot(tvals,epsilonvals) 
plt.show();