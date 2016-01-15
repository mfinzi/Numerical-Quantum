# Author: Marc Finzi
# Last updated: 1/14/2016
# Contactable at mfinzi@hmc.edu

import globalV
import os,sys,time
import FileDialog
import Tkinter as tk
import ttk
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib import patheffects
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from scipy.ndimage import zoom
from scipy.signal import tukey
# Home made imports
from solvers import *
from animator import *
import creationFunctions


class topApp(tk.Tk):

    def __init__(self):

        tk.Tk.__init__(self)

        #self.geometry("650x500+300+300")
        tk.Tk.wm_title(self,"Numerical SE Sim")
        #self.resizable(0,0)
        #tk.Tk.wm_iconbitmap(self,'ih.ico') #Why does this kill the render?

        container = tk.Frame(self, relief="sunken")
        container.pack(side="top", padx=5,pady=5,fill="both",expand=True)
        #for i in np.arange(5):
        #    container.grid_rowconfigure(i,pad=3, weight=1)
        #    container.grid_columnconfigure(i,pad=3, weight=1)

        self.mainPageFrame = MainPage(container,self)
        #self.mainPageFrame.grid(row=0,column=0,sticky="nsew")
        self.mainPageFrame.tkraise()
        

    def quit(self):
        self.destroy()
        sys.exit()

class MainPage(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        
        # The main page frame owns the two figures, the Animator and the Graph communicate
        # through the figures
        self.psiFig = mpl.figure.Figure(figsize=(8,8), dpi=100)#10
        # Initializing the particle into a default config, editable later
        #self.myPart = Particle1d((lambda x: x*x*0.5)(np.linspace(-10,10,2000)),gaussian(0,1,10)(np.linspace(-10,10,2000)),
        #                N=2000,Lx=20.0,dt=.01,SOLVER=EigenSolver1d)
        
        X = np.linspace(-10,10,1024)
        V = creationFunctions.squareBarr(x0=0.1,width=.2,height=100)(X)
        V[np.abs(X)>9.9] = 10000

        self.myPart = Particle1d(V,creationFunctions.gaussian(-5,1,10)(X),
                         N=1024,Lx=20.0,dt=.01,SOLVER=EigenSolver1d)


        # This is the object that will push updates to the figure
        self.anim = Animator(self.myPart,self.psiFig)
        #self.toggleXKCD()
        self.initGrid()

        figContainer = tk.Frame(self, relief="sunken"); figContainer.pack(side=tk.LEFT)
        self.drawSpace = DrawableGraph(self.psiFig,figContainer)
        self.drawSpace.get_tk_widget().pack(side=tk.TOP)#grid(row=1,column=0)#,
                #rowspan=4,padx=5)#,sticky='EN')#pack(side=tk.BOTTOM)

        nonFigContainer = tk.Frame(self, relief="sunken")
        nonFigContainer.pack(side=tk.RIGHT)
        #constantContainer = tk.Frame(nonFigContainer, relief="sunken")

        
        nuclearContainer = tk.Frame(nonFigContainer, relief = "sunken"); nuclearContainer.pack(side = tk.TOP)
        label = tk.Label(nuclearContainer,text="Control Panel", font=("Comic",18)); label.pack()
        quitButton = ttk.Button(nuclearContainer,text="Quit",command=self.controller.quit); quitButton.pack(side=tk.RIGHT)
        resetButton = ttk.Button(nuclearContainer,text='Reset',command = ()); #TODO hook up reset
        resetButton.pack(side = tk.TOP)


        


 
        # Display |psi| or real(psi),im(psi)
        animOptionsContainer = tk.Frame(nonFigContainer, relief = "sunken"); animOptionsContainer.pack()

        subAnimCont = tk.Frame(animOptionsContainer, relief = "sunken"); subAnimCont.pack()
        updateButton = ttk.Button(subAnimCont,text="Redraw",command=self.changeGraph); updateButton.pack(side = tk.LEFT)
        playPauseButton = ttk.Button(subAnimCont,text='Play/Pause',command = globalV.pauseUnpause); playPauseButton.pack(side =tk.RIGHT)

        drawControlsContainer = tk.Frame(animOptionsContainer, relief = "sunken"); drawControlsContainer.pack()
        

        presetContainer = tk.Frame(animOptionsContainer, relief = "sunken"); presetContainer.pack()
        presetLabel = tk.Label(presetContainer,text="Presets", font=("Comic Sans",12)); presetLabel.pack()

        self.presetDictionary = {
                                'Barrier Partial Reflection':("gaussian(-5,1,10)(x)","squareBarr(x0=0.1,width=.2,height=100)(x)"),
                                'Airy functions':("gaussian(-5,1,13)(x)","15*x"),
                                'Harmonic Oscillator':("gaussian(0,2,5)(x)",".5*x**2"),
                                'Abs(x) with Bumps':("(x>-4)*(x<0)*(np.exp(1j*8*x)/np.sqrt(2))*(np.sin(np.pi*x/2))","8*abs(x)+(x<7)*(x>5)*50*(5-x)*(x-7)-(x<-5)*(x>-7)*50*(-5-x)*(-x-7)"),
                                'Coulomb Like':("gaussian(3,1,-12)(x)","-80/(.25*(x-3)**2+.5)-120/(.25*(x+3)**2+.5)")
                                 }

        self.preset = tk.StringVar(); self.preset.set('Barrier Partial Reflection')
        presetsBox = ttk.Combobox(presetContainer,textvariable=self.preset)
        presetsBox['values'] = [key for key in self.presetDictionary]
        presetsBox.pack(side=tk.BOTTOM)
        self.preset.trace('w',self.presetCallback)


        #Todo: Connect up dontInterpBox
        self.dontInterp = tk.BooleanVar()
        dontInterpBox = ttk.Checkbutton(drawControlsContainer,text="Don't Interpo",variable=self.dontInterp);# dontInterpBox.pack(side = tk.RIGHT)

        self.startingK = tk.DoubleVar() # Initial K for psi
        kSlider = ttk.Scale(drawControlsContainer,orient="h",from_=-30,to=30,variable=self.startingK)
        kSlider.pack()

        displayOptionsContainer = tk.Frame(animOptionsContainer, relief = "sunken"); displayOptionsContainer.pack()
        dispTypeButton = ttk.Button(displayOptionsContainer,text='Display Re{'+globalV.uniPsi+'}', command = self.anim.switchDisplayType)
        dispTypeButton.pack(side=tk.LEFT)
        XKCDButton = ttk.Button(displayOptionsContainer, text="XKCD", command=self.toggleXKCD); XKCDButton.pack(side=tk.RIGHT)
        



        # Text inputs for psi and V
        inputContainer = tk.Frame(nonFigContainer, relief="sunken"); inputContainer.pack()
        inputLabel = ttk.Label(inputContainer,text="Direct input",font=("Comic",16)); inputLabel.pack(side=tk.TOP)

        psiContainer = tk.Frame(inputContainer, relief="sunken"); psiContainer.pack()
        psiDLabel = ttk.Label(psiContainer,text=globalV.uniPsi+"(x,0): example exp(-x**2 + 5*1j)"); psiDLabel.pack(side=tk.TOP)
        self.textPsi = tk.StringVar()
        psiInputBox = tk.Entry(psiContainer, textvariable=self.textPsi)
        psiInputBox.pack(side = tk.LEFT)
        self.usePsi = tk.BooleanVar()
        psiCheckBox = ttk.Checkbutton(psiContainer,text="Enable",variable=self.usePsi); psiCheckBox.pack(side = tk.RIGHT)

        vContainer = tk.Frame(inputContainer, relief="sunken"); vContainer.pack()
        vDLabel = ttk.Label(vContainer,text="V(x): example ((x<-5)|(x>5))*100"); vDLabel.pack(side = tk.TOP)
        self.textV = tk.StringVar()
        vxInputBox = tk.Entry(vContainer, textvariable=self.textV)
        vxInputBox.pack(side = tk.LEFT)
        self.useV = tk.BooleanVar()
        vCheckBox = ttk.Checkbutton(vContainer,text="Enable",variable=self.useV); vCheckBox.pack(side = tk.RIGHT)


        #todo add other button functions
        solverContainer = tk.Frame(nonFigContainer, relief = "sunken")
        solverContainer.pack(side=tk.BOTTOM)
        solverTypesContainer = tk.Frame(solverContainer, relief = "sunken")
        solverTypesContainer.pack(side = tk.TOP,expand=True)
        FinDiffButton = ttk.Button(solverTypesContainer, text="Finite Difference", command=(lambda: self.myPart.reInit(SOLVER=EigenSolver1d)))
        FinDiffButton.pack(side = tk.LEFT,fill=tk.BOTH)
        SplitStepButton = ttk.Button(solverTypesContainer, text="Split Step Fourier", command=(lambda: self.myPart.reInit(SOLVER=SplitStepper1d)))
        SplitStepButton.pack(side = tk.RIGHT,fill=tk.BOTH)
        

                # Solver settings
        solverSettingsContainer = tk.Frame(solverContainer, relief = "sunken")
        solverSettingsContainer.pack(side=tk.BOTTOM)
        
        stencilContainer = tk.Frame(solverSettingsContainer); stencilContainer.pack(side = tk.LEFT)
        stencilDescription = ttk.Label(stencilContainer,text="Hamiltonian Stencil"); stencilDescription.pack()
        self.numSTerms = tk.IntVar(); self.numSTerms.set(1)
        stencil3Button = ttk.Radiobutton(stencilContainer,text = "3 Term",variable = self.numSTerms, value = 1); stencil3Button.pack()
        stencil5Button = ttk.Radiobutton(stencilContainer,text = "5 Term",variable = self.numSTerms, value = 2); stencil5Button.pack()
        stencil7Button = ttk.Radiobutton(stencilContainer,text = "7 Term",variable = self.numSTerms, value = 3); stencil7Button.pack()

        #Todo: fix placement of Nscale
        nContainer = tk.Frame(solverSettingsContainer, relief = "sunken"); nContainer.pack(side = tk.RIGHT)
        Ndescription = ttk.Label(nContainer,text="Samples: 1024",); Ndescription.pack(side=tk.TOP)
        self.vectorLength = tk.IntVar()
        self.vectorLength.set(10)
        NScale = ttk.Scale(nContainer,orient='v',from_=3,to=13,variable = self.vectorLength,
            command = lambda x: self.nPointsSliderCallback(Ndescription)); NScale.pack()





        # eigenSlider and coefficients
        altFigSettingsContainer = tk.Frame(figContainer,relief=tk.RAISED,borderwidth=1)
        altFigSettingsContainer.pack(side=tk.BOTTOM,expand=True,fill=tk.BOTH)
        self.eigenNum = tk.IntVar()
        fnum = tk.Label(altFigSettingsContainer)
        fnum.pack(side = tk.LEFT)
        cNum = tk.Label(altFigSettingsContainer)
        cNum.pack(side = tk.RIGHT)
        self.eigenFunctionSlider = ttk.Scale(altFigSettingsContainer,orient="h",from_=0,to=90,variable = self.eigenNum, command = lambda x: self.altGraphSliderCallback((fnum,cNum)) )
        self.eigenFunctionSlider.pack(expand=True,fill=tk.BOTH)

        
    def presetCallback(self,*args,**kwargs):
        self.textPsi.set(self.presetDictionary[self.preset.get()][0])
        self.textV.set(self.presetDictionary[self.preset.get()][1])
        self.usePsi.set(True)
        self.useV.set(True)


    def nPointsSliderCallback(self,label):
        pointExponent = self.vectorLength.get()
        label.config(text = "Samples: %i"%2**pointExponent)

    #def threadedChangeGraph(self):
    #    thread.start_new_thread(self.changeGraph, ())
    def altGraphSliderCallback(self,(fnum,cNum)):
        Ni = self.eigenNum.get()
        fnum.config(text="N_%i"%Ni)
        if self.myPart.solverClass == EigenSolver1d:
            a,b = np.real(self.myPart.solver.basisCoeffs[Ni]),np.imag(self.myPart.solver.basisCoeffs[Ni])
            cNum.config(text="C*C = %.2E"%(a**2 +b**2))
            self.anim.eigenDisplayNum = Ni
            #print np.vdot(self.myPart.solver.basisCoeffs,self.myPart.solver.basisCoeffs)

    def changeGraph(self):
        # Get the drawn curves from the canvas
        dasCurves,updatePsi,updateV = self.drawSpace.extractCurves()

        # downSampling to 50 points with spline interpolation, then rescaling #### NWAS HERE, may have to mess with interpolation settings
        numSamples = 620.
        downSampleRatio = numSamples/float(np.shape(dasCurves)[1])

        rescaleRatio = self.anim.particle.N/float(numSamples)

        pYvY = zoom(dasCurves,(1,downSampleRatio),order = 1)
        pYvY = zoom(pYvY,(1,rescaleRatio),order=1).astype('complex128')


        if self.anim.displayType == 0: # For displaytype and thus drawtype is |psi|^2
            pYvY[0] = np.sqrt(np.abs(pYvY[0]))*np.exp(1j*self.startingK.get()*self.anim.particle.X)
        else:
            pYvY[0] = pYvY[0]*np.exp(1j*self.startingK.get()*self.anim.particle.X)


        oldPsi = self.anim.particle.getPsi()
        oldV = self.anim.particle.getV()
        # We may also need to rescale N if it has changed
        newN = 2**self.vectorLength.get()
        if self.anim.particle.N!= newN:
            scaleRatio = newN/float(self.anim.particle.N)
            oldPsi = zoom(np.real(oldPsi),scaleRatio,order = 1)+1j*zoom(np.imag(oldPsi),scaleRatio,order = 1)
            oldV = zoom(oldV,scaleRatio,order = 1)
            pYvY = zoom(np.real(pYvY),(1,scaleRatio),order=1) +1j*zoom(np.imag(pYvY),(1,scaleRatio),order=1)
 
        newPsi = pYvY[0] if (updatePsi) else oldPsi
        newV   = np.real(pYvY[1]) if (updateV)   else oldV

        # If the direct input boxes are checked, input is taken directly
        X = np.linspace(-self.anim.particle.Lx/2,self.anim.particle.Lx/2,newN) # note that this will inhibit future change of Lx

        if (self.usePsi.get()): 
            lambdaPsi = creationFunctions.getFunctionFromText(self.textPsi.get())
            if lambdaPsi!=None: newPsi = lambdaPsi(X)
        if (self.useV.get()):
            lambdaV = creationFunctions.getFunctionFromText(self.textV.get())
            if lambdaV!=None: newV = lambdaV(X)

        
        # set the new particle settings
        self.anim.particle.reInit(psi = newPsi,Vx = newV,stencilNum = self.numSTerms.get(),N=newN)

        # We need to change the bounds on the altgraph eigenslider
        if self.anim.particle.solverClass == EigenSolver1d:
            newMaxNBasis = self.anim.particle.solver.eigenBasis.shape[1]-1
            self.eigenFunctionSlider.configure(to = newMaxNBasis)
            # For out of bounds safety
            self.anim.eigenDisplayNum = min(self.anim.eigenDisplayNum,newMaxNBasis)            


    def initGrid(self):
        #ttk.Style().theme_use("xpnative")
        self.pack(fill=tk.BOTH,expand=1)

        self.columnconfigure(1, weight=1)
        self.columnconfigure(3, pad=7)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(5, pad=7)



    def toggleXKCD(self):
        plt.xkcd()
        self.psiFig = mpl.figure.Figure(figsize=(8,4), dpi=100)
        self.anim = Animator(self.myPart,self.psiFig)



class Slider(tk.Frame):
    def __init__(self, parent=None ):
        tk.Frame.__init__(self, parent)
        self.number = 0
        self.slide= ttk.Scale(self)
        # self.slide = tk.Scale(self, orient="v", command=self.setValue,
        #                    length=200, sliderlength=20, resolution = .01,
        #                    showvalue=0, tickinterval=2,
        #                    fro=-2, to=2, font=('Arial',9))
        self.text = tk.Label(self, font=('Arial',18))
        self.slide.pack(side=tk.RIGHT, expand=1, fill=tk.X)
        self.text.pack(side=tk.TOP, fill=tk.BOTH)
        #self.unimap = {'4':u'\u2074','5':u'\u2075','6':u'\u2076',
        #               '7':u'\u2077','8':u'\u2078','9':u'\u2079'}

    def setValue(self, val):
        self.number = (10**(int(val)))
        self.text.configure(text='10%s' %0)


# class ButtonDraw():
#     def __init__(self,drawingGraph,master):

#
# Remember to cut window off
#
class DrawableGraph(FigureCanvasTkAgg):
    def __init__(self,figure,master):
        FigureCanvasTkAgg.__init__(self,figure,master)
        self.master = master
        self.cWidth,self.cHeight = figure.get_figwidth()*figure.dpi,\
                      figure.get_figheight()*figure.dpi
        self.b1 = "up"
        self.b2 = "up"
        self.x1old,self.y1old = None,None
        self.x2old,self.y2old = None,None
        self.psiXlist,self.psiYlist = -1*np.ones(self.cWidth),575*np.ones(self.cWidth)
        self.VXlist,self.VYlist = -1*np.ones(self.cWidth),575*np.ones(self.cWidth)
        self.oldPsi = None # stores psYlist and 
        self.oldV = None   # VYlist  for modification
        self.get_tk_widget().bind("<Motion>", self.motion)
        self.get_tk_widget().bind("<ButtonPress-1>", self.b1down) # left click 
        self.get_tk_widget().bind("<ButtonRelease-1>", self.b1up) # for psi

        self.get_tk_widget().bind("<ButtonPress-3>", self.b2down) # right click
        self.get_tk_widget().bind("<ButtonRelease-3>", self.b2up) # for V(x)

        self.get_tk_widget().bind("<ButtonPress-2>",self.loadOldCurves) # So you can extend a curve
        self.get_tk_widget().bind("<space>",lambda event: globalV.pauseUnpause())


    def loadOldCurves(self,event):
        if self.oldPsi != None: self.psiYlist = self.oldPsi
        if self.oldV != None: self.psiYlist = self.oldV


    def extractCurves(self):
        pX,pY,vX,vY = self.psiXlist,self.psiYlist, self.VXlist,self.VYlist
        self.b1 = "up"
        self.b2 = "up"
        self.x1old,self.y1old = None,None
        self.x2old,self.y2old = None,None
        self.psiXlist,self.psiYlist = -1*np.ones(self.cWidth),575*np.ones(self.cWidth)
        self.VXlist,self.VYlist = -1*np.ones(self.cWidth),575*np.ones(self.cWidth)
        self.get_tk_widget().delete("line")
        # Apply tukey window function to psi
        pwY = tukey(620,alpha=.1)*(pY[100:720]-575)/85.
        vwY = vY[100:720]-575
        dasCurves = np.array([pwY,vwY])
        thresholdPsi = np.sum(pX[100:720])>-610
        thresholdV = np.sum(vX[100:720])>-610
        if thresholdPsi: self.oldPsi = pY
        if thresholdV: self.oldV = vY
        return dasCurves, thresholdPsi, thresholdV# Thresholds on activation

    def b1down(self,event):
        if self.inBounds(event.x,event.y):
            self.b1 = "down"           #
            self.get_tk_widget().config(cursor="target") # 
            #globalV.pause()

    def b2down(self,event):
        if self.inBounds(event.x,event.y):
            self.b2 = "down"           #
            self.get_tk_widget().config(cursor="tcross")
            #globalV.pause() # 
            #print event.x,event.y


    def b1up(self, event):
        self.b1 = "up"
        self.x1old = None           # reset the line when you let go of the button
        self.y1old = None
        self.get_tk_widget().config(cursor="arrow")

    def b2up(self, event):
        self.b2 = "up"
        self.x2old = None           # reset the line when you let go of the button
        self.y2old = None
        self.get_tk_widget().config(cursor="arrow")


    def inBounds(self,x,y):
        xGood = (x<self.cWidth) and (x>=0)
        yGood = (y<self.cHeight) and (y>=0)
        return (xGood and yGood)

    def linearInterpRem(self,coords):
        xold,yold,xnew,ynew=coords
        slope = (ynew-yold)/float(xnew-xold)
        return lambda x: self.cHeight - ((x-xold)*slope + yold) #Switch up down

    def motion(self,event):
        if (self.b1 == "down") and self.inBounds(event.x,event.y):
            if self.psiXlist[event.x]==-1:
                color = "black"
            else: color = "red";

            if self.x1old is not None:
                self.get_tk_widget().create_line(self.x1old,self.y1old,event.x,event.y,smooth=True,width=2,fill=color,tag="line")
                if color == "black" and self.x1old!=event.x:
                    coords = self.x1old,self.y1old,event.x,event.y
                    xinRange = np.arange(self.x1old,event.x+1)
                    self.psiXlist[self.x1old:event.x+1] = 1#xinRange
                    self.psiYlist[self.x1old:event.x+1] = self.linearInterpRem(coords)(xinRange)
            self.x1old = event.x
            self.y1old = event.y

        if (self.b2 == "down") and self.inBounds(event.x,event.y):
            if self.VXlist[event.x]==-1:
                color = "blue"
            else: color = "purple";

            if self.x2old is not None:
                self.get_tk_widget().create_line(self.x2old,self.y2old,event.x,event.y,smooth=True,width=2,fill=color, tag="line")
                if color == "blue" and self.x2old != event.x:
                    coords = self.x2old,self.y2old,event.x,event.y
                    xinRange = np.arange(self.x2old,event.x+1)
                    self.VXlist[self.x2old:event.x+1] = 1#xinRange
                    self.VYlist[self.x2old:event.x+1] = self.linearInterpRem(coords)(xinRange)
            self.x2old = event.x
            self.y2old = event.y



if __name__ == "__main__":
    app = topApp()
    app.mainPageFrame.anim.animate()
    app.mainloop()
    
