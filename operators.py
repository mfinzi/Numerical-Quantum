import numpy as np
import scipy.linalg
import scipy as sp
		
class Operator(object):
	""" The most general matrix formulated operator class,
		the other operator classes inherit from this one"""

	def __init__(self,matrix,particle):
		""" Constructor takes in a proper (2darray) matrix"""
		assert matrix.shape[0]==matrix.shape[1] #Must be square
		self.dx = particle.dx
		self.matrix = matrix

	def isHermitian(self):
		"""This is an expensive O(n^2) operation"""
		return self.matrix.H == self.matrix

	def operate(self,vector):
		"""Expensive O(n^2), performs regular matrix mult on a vector"""
		return np.dot(self.matrix,vector)

	def repOperate(self,vector,n):
		"""Calls the operate method n times on vector"""
		assert n>0
		operatedVector = vector
		for i in np.arange(n):
			operatedVector = self.operate(vector)
		return operatedVector

	def expValue(self,vector):
		""" as you would expect, \int_{-\inf}^{\inf} \psi^* A_{op} \psi dx"""
		return np.real(np.vdot(vector,self.operate(vector)))*self.dx

	def uncertainty(self,vector):
		""" as you would expect, \int_{-\inf}^{\inf} \psi^* A_{op} \psi dx"""
		expASqrd = np.vdot(vector,self.operate(self.operate(vector)))*self.dx
		return np.sqrt(expASqrd - self.expValue(vector)**2)

	def computeEigenFunctions(self):
		if self.isHermitian():
			eigenRoutine = sp.linalg.eigh
		else:
			eigenRoutine = sp.linalg.eig

		lambdas,basis = eigenRoutine(self.matrix,overwrite_a=True)
		
		basis/=np.sqrt(particle.dx) # Normalization, see note
		return lambdas, basis

class BandedOp(Operator): #Must also be hermitian

	def __init__(self,bandedM,particle):
		"""Takes in the banded matrix in form [D0,D-1,D-2,...],
		initialization and multiplication is much faster than Operator"""
		assert bandedM.shape[0]%2==1 #Must have odd M
		self.dx = particle.dx
		self.bandedMatrix = bandedM
		

	def isHermitian(self):
		return True

	def computeEigenFunctions(self,jRange=None):
		if jRange: sel='v'
		else: sel ='a' 

		lambdas,basis = sp.linalg.eig_banded(self.bandedMatrix,lower=True, check_finite=False,
					overwrite_a_band=True, select=sel, select_range = jRange)
		basis/=np.sqrt(self.dx) #For normalization
		return lambdas,basis

	def operate(self,vector):
		""" An O(n*m) operation where m is the number of bands and n
		is the size of the vector being operated on"""
		M = self.bandedMatrix.shape[0]
		N = self.bandedMatrix.shape[1]
		#print type(self.bandedMatrix[0][0]),type(self.bandedMatrix[0])
		resultVector = self.bandedMatrix[0]*vector
		for i in np.arange(1,M):
			resultVector += self.shiftDown(np.conj(self.bandedMatrix[i])*vector,i)
			resultVector += self.shiftDown(vector,-i)*self.bandedMatrix[i]
		return resultVector

	def shiftDown(self,vector,n): 
		"""Helper function that fills shifted space with zeros"""
		size = vector.shape[0]
		shifted = np.zeros_like(vector)

		if (n>0):
			shifted[n:] = vector[:size-n] 
		elif (n<0):
			n*=-1
			shifted[:size-n] = vector[n:]
		elif n==0:
			shifted = vector

		return shifted

class XOP(BandedOp):
	def __init__(self,particle):
		self.dx = particle.dx
		self.bandedMatrix = np.reshape(particle.X,(1,particle.N))

class MomentumOp(BandedOp):
	"""Implements the Momentum operator efficiently"""
	# I could implement the stencil table here, but there is no need
	def __init__(self,particle):
		self.dx = particle.dx
		self.bandedMatrix =  np.empty([2,particle.N],dtype=np.complex_)##
		self.bandedMatrix[0] = np.zeros(particle.N)
		self.bandedMatrix[1] = -1j*particle.hbar* np.ones(particle.N)/(2*particle.dx)

class HamiltonianOp(BandedOp):

	STENCIL_TABLE = np.array([-2.0,1.0,-5/2.0,4/3.0,\
					-1/12.0,-49/18.0,3/2.0,-3/20.0,1/90.0])
	SIGMA_CUTOFF = 15

	def __init__(self,particle):
		self.dx = particle.dx
		self.bandedMatrix = self.getDiagonals(particle)
		self.cutOffEnergy = self.expValue(particle.psi)+\
			HamiltonianOp.SIGMA_CUTOFF*self.uncertainty(particle.psi)
	#I could make further optimization by overriding operate with specialized
	# But it may not be necessary
	def getDiagonals(self,P):
		stNum = P.stencilNum
		assert stNum>0, "problem with stencilnum not being positive"
		baseIndex = stNum*(stNum+1)/2 -1 #Stencil table is triangular so we use triangular #'s
		diagArray = np.empty([stNum+1,P.N],dtype=np.float   )
		rho = -P.hbar*P.hbar/(2*P.mass*P.dx*P.dx) #Add some latex equation to explain

		diagArray[0] = P.Vx+rho*HamiltonianOp.STENCIL_TABLE[baseIndex]
		for i in np.arange(1,stNum+1):     
			coeffi = rho*HamiltonianOp.STENCIL_TABLE[i+baseIndex] 
			diagArray[i] = coeffi*np.ones(P.N)		
		return diagArray

	def computeEigenFunctions(self):
		eRange =(-np.inf,self.cutOffEnergy)
		return super(HamiltonianOp,self).computeEigenFunctions(eRange)