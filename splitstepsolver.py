import numpy as np
import scipy.fftpack
import scipy.linalg
import scipy.misc
import scipy as sp
from QM.qmproj.operators import *
import time
import sys
import matplotlib
from matplotlib import pyplot as plt;

class Particle1d:

	def __init__(self,Vx,psi0,N,Lx,dt,SOLVER,consts=[1,1],stencilNum=1,mEAO=1):
		self.stencilNum=stencilNum;
		self.Vx = Vx;
		self.psi = psi0;
		self.psi0 = psi0;
		self.Lx = Lx;
		self.N = N;
		self.dx = float(Lx)/N; ## MAKE SURE Lx IS A FLOAT!!!!!!
		self.X = np.linspace(-Lx/2,Lx/2,N)
		self.mass,self.hbar = consts
		self.dt = dt
		self.t = 0
		self.halt = False
		self.xOp = XOP(self)
		self.momOp = MomentumOp(self)
		self.hamOp = HamiltonianOp(self)
		self.matrExpApproxOrder = mEAO
		
		self.solver = SOLVER(self)
		#try: self.solver = SOLVER(self) #For example self.solver=splitStepper1d(self)
		#except: print "not a valid solver", sys.exc_info()
		
	def reInit(self,**kwargs):
		resize=False
		self.t = 0							# Start at a fresh time

		for key, value in kwargs.iteritems():
			if (key=="consts"):
				self.mass,self.hbar = value
			elif (key=="mEAO"):
				self.matrExpApproxOrder = value
			elif (key!="SOLVER"):			# Catchall for other named keys to
				exec("self."+key +"= value")# update entries with that key

			if (key=="Lx" or key=="N"): 
				resize = True 				# Set resize flag
			if (key=="SOLVER" or key=="stencilNum" or (resize)):
				resolve = True				# Set resolve flag

		if resize:
			self.Lx = float(self.Lx) 		# For safety
			self.dx = self.Lx/self.N; 		# MAKE SURE Lx IS A FLOAT!!!!!!
			self.X = np.linspace(-self.Lx/2,self.Lx/2,self.N)
			self.xOp = XOP(self)			# Since the dimensions changed,
			self.momOp = Momentum(self)		# We need to update the operators
			self.hamOp = Hamiltonian(self)

		if resolve:
			try: self.solver = SOLVER(self) # Recreate the solver
			except: print "not a valid solver",sys.exc_info()

	def step(self):
		self.solver.fullStep(self)

	def getPsi(self):
		return self.psi

	def getPsi_st_psi(self):
		return np.real(np.conj(self.psi)*self.psi)




class BaseSolver(object):
	"""An abstract class detailing the methods to be shared"""
	def fullStep(self,particle):
		return
	def getExpectedEnergy(self,particle):
		return particle.hamOp.expValue(particle.psi)

	def getExpectedMomentum(self,particle):
		return particle.momOp.expValue(particle.psi)

	def getExpectedX(self,particle):
		return particle.xOp.expValue(particle.psi)


class SplitStepper1d(BaseSolver):

	def __init__(self,P):
		self.name="split"
		self.fPsi = None
		# Since <k>=<p>/hbar we can do this
		self.k0 = P.momOp.expValue(P.psi)/P.hbar

		self.dk = 2*np.pi/P.Lx
		self.kx = np.array([np.arange(0,P.N/2),np.arange(-P.N/2,0)]).ravel()*self.dk#+self.k0
		#print self.kx[511],self.kx[512]
	def quarterStpX(self,P):
		# We simply perform the Vx (diagonal) operation in x space
		P.psi *= np.exp(-1j*P.Vx*P.dt/(4*P.hbar))

	def halfStpK(self,P):

		# Recenter Psi in k space, shift fpsi backwards by k0, fpsi(k) -> fpsi(k+k0)
		self.fPsi = np.fft.fft(P.psi*np.exp(-1j*self.k0))*P.dx/np.sqrt(2*np.pi)
		# Then account for that shifting in self.kx by adding in self.k0
		#self.kx -=self.k0
		# Perform the kinetic energy T step in fourier space
		stpFPsi = np.exp(self.kx*self.kx*P.dt*-1j*P.hbar/(4*P.mass))*self.fPsi
		# Then bring kx back to its unshifted array
		
		
		#print self.k0
		#print np.vdot(self.fPsi,self.fPsi)*self.dk
		#self.kx -= self.k0
		#print self.k0
		# Finally we update psi, undoing the kspace shift
		P.psi = np.fft.ifft(stpFPsi)*np.exp(+1j*self.k0)*np.sqrt(2*np.pi)/P.dx

		#self.kx-=self.k0
		# Update k0 while we are in k space 
		self.k0 = self.getExpectedMomentum(P)/P.hbar
		#print self.k0
		
		#self.kx+=self.k0

	def fullStep(self,particle):
		particle.t+=particle.dt
		self.quarterStpX(particle);
		self.halfStpK(particle);
		self.quarterStpX(particle);

	def getExpectedMomentum(self,particle):
		""" Slightly more efficient O(2n) implementation since p is diagonal in kspace"""
		return np.real(np.vdot(self.fPsi,self.fPsi*self.kx))*self.dk


class PropogateSolver1d(BaseSolver):

	def __init__(self,particle):
		self.name="prop"

	def fullStep(self,particle):
		"""Uses the the taylor expansion to the matrix exponential
			for dt << 1, it approximates exp(-itH/hbar)psi =  psi - it/hbar Hpsi +...
			according to the matrExpApproxOrder"""
		particle.t+=particle.dt 			
		temporaryPsi = np.copy(particle.psi)  # We don't want a reference to psi
		for i in np.arange(1,particle.matrExpApproxOrder+1):
			expansionCoeff =(((-1j*particle.dt/particle.hbar)**i)/sp.misc.factorial(i))
			temporaryPsi += expansionCoeff*particle.hamOp.repOperate(particle.psi,i)
		particle.psi = temporaryPsi


class EigenSolver1d(BaseSolver):

	def __init__(self,particle):
		self.name="eigen"
		self.eigenEnergies, self.eigenBasis = self.computeEigenBasis(particle)
		self.basisCoeffs = self.transformToEigenBasis(particle)
		self.reducedEs,self.reducedBasis,self.reducedCoeffs = self.trimBasis()

	def fullStep(self,particle):
		particle.t+=particle.dt
		cT = self.reducedCoeffs*np.exp(-1j*particle.t*self.reducedEs)
		particle.psi = np.dot(self.reducedBasis,cT);

	def computeEigenBasis(self,particle):
		t0 = time.time()
		Es,basis = particle.hamOp.computeEigenFunctions()
		
		print "\n Computation time is %.3fs with %d eigenFunctions used" \
			% ((time.time()-t0),np.size(Es))
		return Es,basis

	def transformToEigenBasis(self,P):
		    #no need for conj here since eigenbasis is real
		basisCoeffs = np.dot(self.eigenBasis.T,P.psi)*P.dx
		return basisCoeffs     #Again, innerProduct = dotProduct * h

	def trimBasis(self):
		arrayMask = np.absolute(self.basisCoeffs)>1e-7
		reducedBasis = self.eigenBasis.T[arrayMask].T
		reducedEs = self.eigenEnergies[arrayMask]#make sure to test this
		reducedCoeffs = self.basisCoeffs[arrayMask]
		print "Trimmed to %d eigenfunctions"%np.size(reducedEs)
		return reducedEs, reducedBasis, reducedCoeffs

	def getExpectedEnergy(self,particle):
		""" Overriding the base class method with a slightly more efficient one"""
		return np.real(np.dot(self.reducedCoeffs*np.conj(self.reducedCoeffs),
			self.reducedEs))
