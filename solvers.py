import time,sys
import numpy as np
from scipy.fftpack import fft, ifft
from scipy.misc import factorial

from operators import *

class Particle1d(object):

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
		self.solverClass = SOLVER
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
			if (key=="SOLVER"): 			# Change the solver type
				self.solverClass = value 
			if (key == "psi"):
				resolve = True
			if (key == "Vx"):
				resize = True
				resolve = True

		if resize:

			self.Lx = float(self.Lx) 		# For safety
			self.dx = self.Lx/self.N; 		# MAKE SURE Lx IS A FLOAT!!!!!!
			self.X = np.linspace(-self.Lx/2,self.Lx/2,self.N);
			#self.Vx[np.abs(self.X)>9.9]=10000     # Add hoc modification to prevent wraparound
			self.xOp = XOP(self)			# Since the dimensions changed,
			self.momOp = MomentumOp(self)		# We need to update the operators
			self.hamOp = HamiltonianOp(self)
			print "resize"

		if resolve:
			try: self.solver = self.solverClass(self) # Recreate the solver
			except: print "not a valid solver",sys.exc_info()

	def step(self):
		self.solver.fullStep(self)

	def getPsi(self):
		return self.psi

	def getV(self):
		return self.Vx

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
	""" Implements the Split Step Fourier method, which is well explained in Jake VanderPlas's blog
		at https://jakevdp.github.io/blog/2012/09/05/quantum-python/ """
	def __init__(self,P):
		self.name="split"
		self.fPsi = None
		# Since <k>=<p>/hbar we can do this
		self.k0 = P.momOp.expValue(P.psi)/P.hbar

		self.dk = 2*np.pi/P.Lx
		self.kx = np.array([np.arange(0,P.N/2),np.arange(-P.N/2,0)]).ravel()*self.dk#+self.k0
		self.stepsPerStep = 20
		self.subDt = P.dt/self.stepsPerStep
		#Precompute the halfStepX exponential
		self.halfXOperation = np.exp(P.Vx*(-1j*self.subDt/(2*P.hbar)))
		#Precompute the quarterStepK exponential
		self.quarterKOperation = np.exp((self.kx**2)*(self.subDt*-1j*P.hbar/(8*P.mass)))
	def halfStpX(self,P):
		# We simply perform the Vx (diagonal) operation in x space
		P.psi *= self.halfXOperation

	def quarterStpK(self,P):

		# Recenter Psi in k space, shift fpsi backwards by k0, fpsi(k) -> fpsi(k+k0)
		self.fPsi = fft(P.psi*np.exp(-1j*self.k0))/np.sqrt(2*np.pi)
		# Then account for that shifting in self.kx by adding in self.k0
		#self.kx -=self.k0
		# Perform the kinetic energy T step in fourier space
		self.fPsi *= self.quarterKOperation
		# Then bring kx back to its unshifted array
		
		
		#print self.k0
		#print np.vdot(self.fPsi,self.fPsi)*self.dk
		#self.kx -= self.k0
		#print self.k0
		# Finally we update psi, undoing the kspace shift
		P.psi = ifft(self.fPsi)*np.exp(+1j*self.k0)*np.sqrt(2*np.pi)

		#self.kx-=self.k0
		# Update k0 while we are in k space 
		self.k0 = self.getExpectedMomentum(P)/P.hbar
		#print self.k0
		
		#self.kx+=self.k0

	def fullStep(self,particle):
		for i in np.arange(self.stepsPerStep):
			particle.t+=self.subDt
			self.quarterStpK(particle);
			self.halfStpX(particle);
			self.quarterStpK(particle);

	def getExpectedMomentum(self,particle):
		""" Slightly more efficient O(2n) implementation since p is diagonal in kspace"""
		return np.real(np.vdot(self.fPsi,self.fPsi*self.kx))*self.dk


class PropogateSolver1d(BaseSolver):
	""" There is some numerical instability, need to investigate further"""
	def __init__(self,particle):
		self.name="prop"

	def fullStep(self,particle):
		"""Uses the the taylor expansion to the matrix exponential
			for dt << 1, it approximates exp(-itH/hbar)psi =  psi - it/hbar Hpsi +...
			according to the matrExpApproxOrder"""
		particle.t+=particle.dt 			
		temporaryPsi = np.copy(particle.psi)  # We don't want a reference to psi
		for i in np.arange(1,particle.matrExpApproxOrder+1):
			expansionCoeff =(((-1j*particle.dt/particle.hbar)**i)/factorial(i))
			temporaryPsi += expansionCoeff*particle.hamOp.repOperate(particle.psi,i)
		particle.psi = temporaryPsi


class EigenSolver1d(BaseSolver):
	""" Performs the decomposition of the wavefunction into the basis of the first N
		eigenfunctions. Essentially performs the diagonalization of the Hamiltonian
		operator, discretized by a finite difference approximation in the Operators file.
		Time evolution is computed from the matrix exponential of the Hamiltonian matrix"""

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
