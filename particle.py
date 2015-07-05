#!/usr/bin/python
import numpy as np;
import scipy as sp;
import matplotlib as mpl
mpl.use("TKAgg")
from matplotlib import pyplot as plt;
from matplotlib import animation;
import time, sys, getopt
import scipy.sparse.linalg

#   Author: Marc Finzi

#   Instructions: This is the source code, either use the compiled version
#   or compile on command line by running python and then running
#   import py_compile;py_compile.compile('particle.py') to generate the pyc file
#   if you already have the pyc file then just run{ python particle.pyc [p0-p5] }
#   you can save the video with an additional -s [filenameHere.mp4] but you must
#   have ffmpeg encoder and the H.264 codec. Results may vary based on OS

class Particle:
    def __init__(self,potential,psi0,h,xMax,firstN=10,xMin='-xMax'):
        self.Vx = potential;
        
        if xMin=='-xMax': xMin = -xMax #Gets around asymmetric bounds
        
        self.eigenEnergies, self.eigenBasis, self.X \
            = self.computeEigenBasis(potential,h,xMax,firstN,xMin);
            
        self.psi = psi0
        
        self.basisCoeffs = self.transformToEigenBasis(h)
        self.h = h
        
        
    def computeEigenBasis(self,Vx,h,xMax,firstN,xMin):
        """ This function finds all of the energy eigenvalues and the 
        corresponding eigenfunctions of schrodinger's equation with potential 
        Vx using the finite difference eigenvalue method. Takes in the x step 
        h, the maximum x value xMax and the function Vx."""  
        t0 = time.time()
        
        xArr = np.arange(xMin,xMax,h)
        numPoints = int(round((xMax-xMin)/h))
        mainDiag = 2*h*h*Vx(xArr)+2;
        lowerDiag = np.zeros(numPoints-1)-1;
        upperDiag = lowerDiag;
        cMatrix = np.diag(mainDiag)+np.diag(lowerDiag,-1)+np.diag(upperDiag,1)
        
        if (firstN>.2*numPoints):
            EsAll,basisAll = np.linalg.eigh(cMatrix); #Use regular method
            Es,basis = EsAll[:firstN],basisAll[:,:firstN]
            
        else: # If we dont want very many Basis functions, 
                         #we can use the sparse approach
            Es,basis = sp.sparse.linalg.eigsh(cMatrix,firstN,sigma=0,which='LM')
            
        basis/=np.sqrt(h) #Although the eigenfunctions are orthonormal vectors, they
                # need to be normalized as functions, dot product misses an h
        print "\n Computation time is %.2fs" % (time.time()-t0)
        return Es/(2*h*h),basis,xArr
    

    def transformToEigenBasis(self,h):
        basisCoeffs = np.dot(self.eigenBasis.transpose().conj(),self.psi)*h
        return basisCoeffs;     #Again, innerProduct = dotProduct * h
    
    def removeZerosFromBasis(self):
    	zeroS = np.where()

    def updatePsi(self,Psi):
        self.psi = Psi;   
        
    def computePsiT(self,t):
        cT = self.basisCoeffs*np.exp(-1j*t*self.eigenEnergies)
        psiT = np.dot(self.eigenBasis,cT)
        return psiT

    def computePsiStarPsi(self,t): # Make sure eigenEnergies are in the same 
        psiT = self.computePsiT(t)  # order as the eigenFunctions
        return np.real(psiT*np.conj(psiT)) # There may be roundoff errors, so take Real
        
    def getExpectedEnergy(self):
        return np.real(np.dot(self.basisCoeffs*np.conj(self.basisCoeffs),self.eigenEnergies))

        
class Animation:
    def __init__(self,particle):
        self.fig = plt.figure()
        
        peakPsi = max([-1*np.min(np.absolute(particle.psi)),np.max(np.absolute(particle.psi))])**2
        # Graph is scaled based on the initial height of psi
        
        self.axe = plt.axes(xlim = (particle.X[0],particle.X[-1]),
                             ylim = (-3*peakPsi,3*peakPsi)) 
        self.line, = self.axe.plot([], [], lw=1) # '-o' optional
        self.axe.add_line(mpl.lines.Line2D(particle.X,
                                    .008*particle.Vx(particle.X), color = 'k'))
        self.axe.grid()
        self.particle = particle
        self.txtBox = self.axe.text(.02, 0.90, '', transform=self.axe.transAxes)
        
    def inits(self):
        
        self.txtBox.set_text('')
        self.line.set_data([], [])
        return self.line, self.txtBox
        
    def update_line(self,time):
        x = self.particle.X
        y = self.particle.computePsiStarPsi(time*.01)
        self.line.set_data(x,y)
        self.txtBox.set_text('<E> = %.3f' %self.particle.getExpectedEnergy())
        return self.line, self.txtBox
    
    def animate(self):
        a = animation.FuncAnimation(self.fig, self.update_line,
                    init_func=self.inits, frames=2000, interval=25, blit=True,repeat = False)
        return a
        
    def h264save(self, anim, name):
        t0 = time.time()
        anim.save(str(name),writer="ffmpeg_file", fps=30,dpi=96,extra_args = ['-vcodec', 'libx264'])
        print "Encode time is %.2fs"%(time.time()-t0)
        
def checkOrthoNormal(P):
    sizeBasis = np.shape(P.eigenBasis)[1]
    identity = np.diag(np.ones(sizeBasis)) #Again the h factor is from innerproduct
    diff = np.dot(P.eigenBasis.transpose(),P.eigenBasis)*P.h - identity
    return np.all(diff*diff)<1E-20 # These are floats, so we will never get exactly 0
            
def squareBarr(x0=0,width=10,height=1E10):
    """Barrier if height is +, Well if -. Centered on x0"""
    return lambda x: (x>(x0-width/2.0))*(x<(x0+width/2.0))*height

def lennardJonesPot(rm=5,depth =400):
    return lambda r: depth*((rm/r)**12 - 2*(rm/r)**6)
   
def gaussian(x0=0,width=1,vel=0):
    """Returns a normalized gaussian function in x,
        please make sure width is nonzero """
    norm = np.sqrt(np.sqrt(np.pi)*float(width))
    return lambda x: np.exp(-.5*((x-x0)/width)**2+vel*1j*x)/norm
    



## THese two parts should really be in a different file
## but for convenience they both in the same file



# Will probably run a bit slow on a laptop

P0 = "Particle(lambda x: np.abs(x)+(x<15)*(x>10)*11*(10-x)*(x-15)-(x<-10)*(x>-15)*11*(-10-x)*(-x-15) , (lambda x2: (x2>-5)*(x2<0)*((np.exp(1j*12*x2)/np.sqrt(2))*(np.sin(np.pi*x2/5))))(np.arange(-20,20,.01)), .01,20,200)"
P0Text = "Abs(x) potential with 2 inverted quadratic barriers"

P1 = "Particle(lambda x: -80/(.25*(x-3)**2+.5)-120/(.25*(x+3)**2+.5), gaussian(3,1,-12)(np.arange(-15,15,.01)),.01,15,3000)"
P1Text = "Asymmetric, double coulomb ish, this one can be problematic,\
   most of the energy eigenvalues are negative and we can only\
   find the first 3000 with this method"

P2 = "Particle(squareBarr(12.5,2.5,80), gaussian(5,vel=14)(np.arange(-0,20,.01)), .01,20,200,0)"
P2Text = "Partial reflection off of a square barrier"

P3 = "Particle(lennardJonesPot(8), gaussian(8,.25,-5)(np.arange(4,20,.005)), .005,20,1000,4)"
P3Text = "Atom trapped in lennard jones"

P4 = "Particle(lambda x: .5*x**2  , gaussian(0,2,5)(np.arange(-15,15,.01)), .01,15,200)"
P4Text = "Your standard sho potential. Interesting behavior when the width of the gaussian is 1, then\
    amplitiude appears not to change. Probably related to the fact that [gaussian(0,1,0)]\
    is an eigenfunction."
    
P5 = "Particle(lambda x: x*10, gaussian(-5,1,13)(np.arange(-15,15,.01)), .01,15,200)"
P5Text =  "Linear potential, looking at P.eigenBasis[:,5] or P.eigenBasis[:,100]... you see the Airy functions"
 


def main(args):

    if len(args)==0: 
        print "Usage: particle.pyc P<0-5> [-s] [name]\n or \nparticle.py manual [-s] [name]"
        sys.exit(0)

    if args[0] == "manual":
        funcArgs = eval(raw_input("type expression for [<potential>, <psi0>, <h>, <xMax>, <firstN>, <xMin>]"))
        P = Particle(*funcArgs)
    else:
        Pn = args[0].upper()
        try: theP = eval(Pn)
        except NameError: print "not an option, choose P0-P5"; sys.exit(0);
        theText = eval(Pn+"Text")
        print theText
        exec("P="+theP) #Executes the chosen particle line

    plt.ioff()
    A = Animation(P);
    anim = A.animate();

    if ("-s" in args): 
        try: saveName = args[2] #Try to get the filename
        except IndexError: saveName = "QM_SIM.mp4" #Default name
        A.h264save(anim,saveName) #Save the animation if -s is found
    else: plt.show()              #and in that case don't display it

    return 0

if __name__ == "__main__":
    main(sys.argv[1:])

