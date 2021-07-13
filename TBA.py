""" class for Thermodynamic Bethe Ansatz calculations """

#libraries
import numpy as np
import parameters

epsilon = 0.0001	#minimal difference between Fermi points
					#to be considered different
N = 50				#number of point used in each interval
					#for computation of discrete integrals

#-------------------------------------------------------------------------
#	Lieb-Liniger kernel
#-------------------------------------------------------------------------

Phi = lambda k1, k2: 2*np.arctan((k1-k2)/parameters.c)
varphi = lambda k1, k2: 2*parameters.c/((k1-k2)**2 + parameters.c**2)

#-------------------------------------------------------------------------
#	definition of Split Fermi Sea class
#-------------------------------------------------------------------------

class SFS(object):
	#construct object 'split Fermi sea' associated
	#with set of Fermi points 'theta'
    def __init__(self, theta):
        self.theta = theta
        self.sigma = []
        #count components and 'special points'
        comp=0; sp=0;
        for j in range(len(theta)/2):
            if theta[2*j+1] > theta[2*j]+epsilon:
                comp+=1
            else:
                sp+=1
            self.sigma.append(-1)
            self.sigma.append(1)
		#properties associated with the object 'Split Fermi Sea'
        self.lam_discr = np.zeros((comp*N+2*sp))	#discrete space of rapidities
        self.varphi_mat = np.zeros((self.lam_discr.shape[0],self.lam_discr.shape[0]))	#Lieb-Liniger kernel matrix
        self.dresser = np.zeros((self.lam_discr.shape[0],self.lam_discr.shape[0]))		#dressing matrix
        self.dlam = np.zeros((self.lam_discr.shape[0]))		#measure on rapidity space (for integration)
        self.indx=[]
        self.veff=np.zeros((len(self.theta)))		#effective velocity at all Fermi points
		#discretized occupation functio
        count=0;
        for j in range(len(theta)/2):
            if self.theta[2*j+1] > self.theta[2*j]+epsilon:
                self.lam_discr[count:count+N] = np.linspace(self.theta[2*j],self.theta[2*j+1],N)
                self.dlam[count+1:count+N-1] = self.lam_discr[count+1:count+N-1]-self.lam_discr[count:count+N-2]
                self.dlam[count] = self.dlam[count+1]/2.
                self.dlam[count+N-1] = self.dlam[count]
                self.indx.append(count)
                self.indx.append(count+N-1)
                count += N
            else:
                self.lam_discr[count:count+2] = np.array([self.theta[2*j]-0.5*epsilon,self.theta[2*j]+0.5*epsilon])
                self.dlam[count:count+2] = epsilon/2.
                self.indx.append(count)
                self.indx.append(count+1)
                count += 2
        for i in range(self.lam_discr.shape[0]):
            for j in range(self.lam_discr.shape[0]):
                self.varphi_mat[i,j] = varphi( self.lam_discr[i], self.lam_discr[j] )
        #calculate dressing matrix
        self.dresser = np.eye(self.lam_discr.shape[0]) - 1/(2*np.pi) * self.varphi_mat.dot(np.diag(self.dlam))
    def dress(self, f):
        h = np.array( [f(self.lam_discr[i]) for i in range(self.lam_discr.shape[0])] )
        return np.linalg.solve(self.dresser, h)
    def charge_density(self, f):
        fdr = self.dress(f)
        return fdr.dot(self.dlam)/(2*np.pi)		
    def calculate_veff(self):
		lamdr = self.dress(lambda u: u)
		onedr = self.dress(lambda u: 1.)
		for j in range(len(self.indx)):
			self.veff[j] = lamdr[self.indx[j]]/onedr[self.indx[j]]