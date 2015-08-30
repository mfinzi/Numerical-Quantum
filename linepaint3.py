import Tkinter as tk
import ttk
import tkMessageBox
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from NumericalQuantum.splitstepsolver import *
from NumericalQuantum.animator import *
from NumericalQuantum.creationfunctions import *

class topApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        #self.geometry("650x500+300+300")
        tk.Tk.wm_title(self,"Numerical SE Sim")

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

class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # The main page frame owns the two figures, the Animator and the Graph communicate
        # through the figures
        self.psiFig = mpl.figure.Figure(figsize=(8,4), dpi=100)
        self.alternateFig = mpl.figure.Figure(figsize=(8,4), dpi=100)

        # Initializing the particle into a default config, editable later
        myPart = Particle1d((lambda x: x*x*0.5)(np.linspace(-10,10,100)),gaussian(0,1,10)(np.linspace(-10,10,100)),
                        N=100,Lx=20.0,dt=.01,SOLVER=EigenSolver1d)

        # This is the object that will push updates to the figure
        anim = Animator(myPart,self.psiFig,self.alternateFig)

        self.initGrid()

        label = tk.Label(self,text="THis is it!", font=("Calibri",12))
        quitButton =ttk.Button(self,text="Quit",command=self.controller.quit)
        #bButton =ttk.Button(self,text="Redraw",command=anim.threadedAnimate())
        

        inputContainer = tk.Frame(self, relief="sunken")
        textPsi = tk.StringVar()
        psiInputBox = tk.Entry(inputContainer, textvariable=textPsi)
        textVx = tk.StringVar()
        vxInputBox = tk.Entry(inputContainer, textvariable=textVx)



        constantContainer = tk.Frame(self, relief="sunken")
        mass = tk.DoubleVar()
        massSlider = Slider(self)#tk.Scale(constantContainer,orient="v",from_=.01,to=100)#,label="Mass")


        solverContainer = tk.Frame(self, relief = "sunken")
       # FinDiffButton = ttk.Button(solverContainer, text="Finite Diff", command=)



        solverSettingsContainer = tk.Frame(self, relief = "sunken")





        drawSpace = DrawableGraph(self.psiFig,self)
        secondaryBox = FigureCanvasTkAgg(self.alternateFig,self)



        #Packing area
        constantContainer.pack()
        massSlider.pack()
        label.pack()
        quitButton.pack()
        inputContainer.pack()
        psiInputBox.pack()
        vxInputBox.pack()
        drawSpace.get_tk_widget().pack()#grid(row=1,column=0)#,
                #rowspan=4,padx=5)#,sticky='EN')#pack(side=tk.BOTTOM)
        #anim.threadedAnimate()

    def initGrid(self):
        #ttk.Style().theme_use("xpnative")
        self.pack(fill=tk.BOTH,expand=1)

        self.columnconfigure(1, weight=1)
        self.columnconfigure(3, pad=7)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(5, pad=7)

    def getFunctionFromText(self,textVariable):
        # I should put this in a different file with #from numpy import *
        arrayString ='(lambda x:'+textVariable+')(self.X)'

        disallowedStrings = [";","()","sys"]
        for antiString in disallowedStrings:
            if antiString in arrayString:
                tkMessageBox.showerror("Sanitizer",
                    "No funny business with the text input")
                return None
        try:
            function = eval(arrayString)
        except SyntaxError:
            tkMessageBox.showwarning("Syntax Error",
                "Base python expressions and numpy functions only,\n"\
                "you don\'t need to worry about normalization for psi. \n"\
                "Example gaussian psi:     50*exp(-(x-3)**2/3 + 5*1j)\n"\
                "Example Small Square well:    ((x<-5)|(x>5))*100\n"\
                "alternate way for making the well, (it goes down this time):\n"\
                " (x>-5)*(x<5)*-100")
            return None
        return function


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


class DrawableGraph(FigureCanvasTkAgg):
    def __init__(self,figure,master):
        FigureCanvasTkAgg.__init__(self,figure,master)
        self.master = master
        self.cWidth,self.cHeight = figure.get_figwidth()*figure.dpi,\
                      figure.get_figheight()*figure.dpi
        self.b1 = "up"
        self.xold,self.yold = None,None
        self.xlist,self.ylist = -1*np.ones(self.cWidth),np.zeros(self.cWidth)

        self.get_tk_widget().bind("<Motion>", self.motion)
        self.get_tk_widget().bind("<ButtonPress-1>", self.b1down)
        self.get_tk_widget().bind("<ButtonRelease-1>", self.b1up)

        #self.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        

    def b1down(self,event):
        if self.inBounds(event.x,event.y):
            self.b1 = "down"           #
            self.get_tk_widget().config(cursor="target") # 
            print event.x,event.y

    def b1up(self, event):
        self.b1 = "up"
        self.xold = None           # reset the line when you let go of the button
        self.yold = None
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
        if self.b1 == "down" and self.inBounds(event.x,event.y):
            if self.xlist[event.x]==-1:
                color = "black"
            else: color = "red";

            if self.xold is not None and self.xold!=event.x:
                self.get_tk_widget().create_line(self.xold,self.yold,event.x,event.y,smooth=True,width=2,fill=color)
                if color == "black":
                    coords = self.xold,self.yold,event.x,event.y
                    xinRange = np.arange(self.xold,event.x+1)
                    self.xlist[self.xold:event.x+1] = xinRange
                    self.ylist[self.xold:event.x+1] = self.linearInterpRem(coords)(xinRange)
            self.xold = event.x
            self.yold = event.y

if __name__ == "__main__":
    app = topApp()
    
    app.mainloop()
    #x = app.mainPageFrame.boss.xlist
    #y = app.mainPageFrame.boss.ylist
    #plt.plot(y)
    #plt.show()
