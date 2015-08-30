import numpy as np;
import matplotlib as mpl
from matplotlib import animation;
mpl.use("TkAgg")
import time
from QM.qmproj import splitstepsolver
import thread
from matplotlib import pyplot as plt;


class Animator:
    def __init__(self,particle,mainFig,altFig):
    	self.particle = particle

    	# Setup the figures
        self.psiFig = mainFig
        self.alternateFig = altFig
        self.psiAxe = self.psiFig.add_subplot(111)
        self.alternateAxe = self.alternateFig.add_subplot(111)

        # Graph is scaled based on the initial height of psi, and the x bounds
        peakPsi = max([-1*np.min(np.absolute(particle.psi)),np.max(np.absolute(particle.psi))])**2
        self.psiAxe.set_xlim(particle.X[0],particle.X[-1])
        self.psiAxe.set_ylim(-3*peakPsi,3*peakPsi)
        
        # Scaling of the alternate graph depends on what type of solver it is using
        # TODO add the scaling for the alternate graph

        # Setup the empty line on the psigraph to hold psi*psi, and also the Vx on that graph
        self.psiLine, = self.psiAxe.plot([], [], lw=1) # '-o' optional
        self.psiVxLine = self.psiAxe.plot(particle.X, .008*particle.Vx,color = 'k')
        self.psiAxe.grid()
        


        self.anim=None
        #self.txtBox = self.axe.text(.02, 0.90, '', transform=self.axe.transAxes)   
    def inits(self):
        
        #self.txtBox.set_text('')
        self.psiLine.set_data([], [])
        plt.plot(0) #Need to figure out why I need this line
        return self.psiLine,
        
    def update_line(self,t):
        x = self.particle.X
        #y = self.particle.computePsiStarPsi(time*.01)
        
        self.particle.solver.fullStep(self.particle)
        
        #print "this step is %.7f" %(time.time()-self.ass)
        y= np.real(self.particle.getPsi_st_psi())
        #self.ass = time.time()
        self.psiLine.set_data(x,y)
        #print "this step is %.7f" %(time.time()-t0)
        #self.txtBox.set_text('<E> = %.3f' %self.particle.getExpectedEnergy())
        return self.psiLine,#, self.txtBox
    
    def threadedAnimate(self):
        thread.start_new_thread(self.animate, ())

    def animate(self):
        self.anim = animation.FuncAnimation(self.fig, self.update_line,
                    init_func=self.inits, frames=100, interval=16, blit=True,repeat = False)
        

    def h264save(self, anim, name):
        t0 = time.time()
        anim.save(name,writer="ffmpeg_file", fps=30,dpi=96,extra_args = ['-vcodec', 'libx264'])
        print "Encode time is %.2fs"%(time.time()-t0)