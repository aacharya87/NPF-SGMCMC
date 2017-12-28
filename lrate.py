import os
import sys
import numpy as np
from pylab import *
from math import *
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


aval = 1.0;
bval = aval;
cval = aval;
tvals = np.arange(100);

epsilonvals = pow(aval*(1+tvals/bval),-cval);

plt.plot(tvals,epsilonvals) 
plt.show();