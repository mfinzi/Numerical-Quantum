import numpy as np
import time

# File Defines some convenient functions to be used by the raw input panel, and presets

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

import tkMessageBox
from numpy import *
#from numpy import linspace
def getFunctionFromText(textVariable):
    # I should put this in a different file with #from numpy import *
    arrayString ='(lambda x:'+textVariable+')'

    disallowedStrings = [";","()","sys"]
    for antiString in disallowedStrings:
        if antiString in arrayString:
            tkMessageBox.showerror("Sanitizer",
                "No funny business with the text input")
            return None
    function = eval(arrayString)
    # test it out
    X = linspace(-10,10,100)
    try: function(X)
    except:
        tkMessageBox.showwarning("Syntax Error",
            "Base python expressions and numpy functions only,\n"\
            "Example gaussian psi:     50*exp(-(x-3)**2/3 + 5*1j)\n"\
            " Or you could use the gaussian function: gaussian(x0=3,width=sqrt(3),vel=5)\n"
            "Example Small Square well:    ((x<-5)|(x>5))*100\n"\
            " Or you could use the function: squareBarr(x0=0,width=10,height=-100)\n"
            "alternate way for making the well, (it goes down this time):\n"\
            " (x>-5)*(x<5)*-100")
        return None
    return function