import numpy as np
import time


def timer(function,n):
	elapsed = 0
	for i in np.arange(n):
		t0 = time.time()
		function()
		dt= time.time()-t0
		elapsed+=dt
	return elapsed/n

def gaussian(x0=0,width=1,vel=0):
    """Returns a normalized gaussian function in x,
        please make sure width is nonzero """
    norm = np.sqrt(np.sqrt(np.pi)*float(width))
    return lambda x: np.exp(-.5*((x-x0)/width)**2+vel*1j*x)/norm

def squareBarr(x0=0,width=10,height=1E10):
    """Barrier if height is +, Well if -. Centered on x0"""
    return lambda x: (x>(x0-width/2.0))*(x<(x0+width/2.0))*height
