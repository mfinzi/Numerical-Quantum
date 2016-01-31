import time #, thread
import numpy as np;

from matplotlib import animation;
from matplotlib import pyplot as plt
from matplotlib import gridspec

import globalV
import solvers



class Animator(object):
    def __init__(self,particle,mainFig):

    	self.particle = particle

    	# Setup the figures
        
        self.psiFig = mainFig
        grid = gridspec.GridSpec(2,1,height_ratios = [4,3])
        self.psiAxe = self.psiFig.add_subplot(grid[0])#211)
        self.altAxe = self.psiFig.add_subplot(grid[1])#212)

        # Graph is scaled based on the initial height of psi, and the x bounds
        peakPsi = max([-1*np.min(np.absolute(particle.psi)),np.max(np.absolute(particle.psi))])**2
        self.altLims = ((particle.X[0],particle.X[-1],-3*peakPsi,3*peakPsi),(-30,30,0,100))
        self.psiAxe.set_xlim(particle.X[0],particle.X[-1])
        self.psiAxe.set_ylim(-3*peakPsi,3*peakPsi)
        self.altAxe.set_xlim(particle.X[0],particle.X[-1])
        self.altAxe.set_ylim(-3*peakPsi,3*peakPsi)
        # Scaling of the alternate graph depends on what type of solver it is using
        #self.altAxe.set

        # Setup the empty line on the psigraph to hold psi*psi, and also the Vx on that graph
        self.psiLine, = self.psiAxe.plot([], [], lw=1) # '-o' optional
        self.psiVxLine, = self.psiAxe.plot(particle.X, .01*particle.Vx,color = 'k')
        self.psiAxe.grid()
        #
        self.altLine, = self.altAxe.plot([],[],lw=2)
        self.altAxe.grid()

        self.anim=None
        #self.txtBox = self.axe.text(.02, 0.90, '', transform=self.axe.transAxes)
        self.displayType = 0
        
        self.eigenDisplayNum = 0
        self.solverOld = "eigen"

    def inits(self):

        #self.txtBox.set_text('')
        self.psiLine.set_data([], [])
        self.psiVxLine.set_data([], [])
        self.altLine.set_data([], [])
        plt.plot(0) #Need to figure out why I need this line
        return self.psiLine,self.psiVxLine,self.altLine,
        
    def update_line(self,t):

        x = self.particle.X
        #y = self.particle.computePsiStarPsi(time*.01)
        
        if not globalV.paused:
            self.particle.solver.fullStep(self.particle)
            #print "this step is %.7f" %(time.time()-self.ass)
        if self.displayType == 0:
                y= np.real(self.particle.getPsi_st_psi())
        else: y = np.real(self.particle.getPsi())
    
        self.psiLine.set_data(x,y)
        self.psiVxLine.set_data(x,.01*self.particle.Vx)
        #print "this step is %.7f" %(time.time()-t0)
        #self.txtBox.set_text('<E> = %.3f' %self.particle.getExpectedEnergy())


        if self.particle.solver.name=="eigen":
            if self.solverOld != "eigen":
                self.altAxe.set_xlim(self.altLims[0][0],self.altLims[0][1])
                self.altAxe.set_ylim(self.altLims[0][2],self.altLims[0][3])
            self.altLine.set_data(x,self.particle.solver.eigenBasis[:,self.eigenDisplayNum])

        if self.particle.solver.name== "split" and not globalV.paused:
            k = self.particle.solver.kx
            fPsi = self.particle.solver.fPsi
            if self.displayType==0:
                y=np.real(fPsi*np.conj(fPsi))
            if self.displayType==1:
                y=np.real(fPsi)
            #normalize y
            #y = y/np.sum(y)
            if self.solverOld != "split":
                self.altAxe.set_xlim(-1*k[np.shape(k)[0]/2]/8,k[np.shape(k)[0]/2]/8)
                self.altAxe.set_ylim(2*np.min(y),2*np.max(y))

            self.altLine.set_data(k,y)

        self.solverOld=self.particle.solver.name
        return self.psiLine,self.psiVxLine,self.altLine,#, self.txtBox
    
    #def threadedAnimate(self):
    #    thread.start_new_thread(self.animate, ())

    def animate(self):
        self.anim = animation.FuncAnimation(self.psiFig, self.update_line, #true
                    init_func=self.inits, frames=1000000, interval=33.3, blit=True,repeat = False)

        

    def h264save(self, anim, name):
        t0 = time.time()
        anim.save(name,writer="ffmpeg_file", fps=30,dpi=96,extra_args = ['-vcodec', 'libx264'])
        print "Encode time is %.2fs"%(time.time()-t0)


    def switchDisplayType(self):
        self.displayType^=1