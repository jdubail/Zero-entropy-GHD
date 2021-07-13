""" file with parameters for zero-entropy GHD simulation """

#libraries
import numpy as np

###########################################################
#
#	parameters
#
###########################################################

kB=1.38065e-23		#Boltzmann constant		(um^2.ms^-2.kg.uK^-1)
hbar = 1.05457e-25	#Planck's constant		(um^2.kg.ms^-1)


mu0= 0.014602 * kB	#chemical potential		(uK * kB)
m=1.44e-25			#mass    				(kg)

c=4.966				#contact interaction (um^-1)

print 'c = ', c, 'um^-1'


x0=35.
l=70.

w=55.17 			#gaussian beam width (um)
omegaanti=41.6e-3 	#anti trap (kHz)

U0=1.504e-24   		#initial trap amplitude 	(uJ/kB=uK)
V0 = lambda z: ((U0*(1.-np.exp(-2*(z/w)**2))- 0.5*m*omegaanti**2*z**2) - mu0) * m/hbar**2

Uf=1.464e-22		#final trap amplitude 		(uJ/kB=uK)
V = lambda z: (Uf*(1.-np.exp(-2*(z/w)**2))- 0.5*m*omegaanti**2*z**2) * m/hbar**2	#potential V(x) for evolution


xspace = np.linspace(-60.,60.,2000)			#discretization of space (points where density is calculated)

omega = np.sqrt((4*Uf/w**2)/m)
print 'omega =', omega
dt = 0.0002 * 2*np.pi/omega		#time step
Delta_t = 0.01 * 2*np.pi/omega	#save files every Delta_t

print 'dt =', dt
print 'Delta_t =', Delta_t
tmax = 10*2*np.pi/omega			#final time

###########################################################
#
#	turns into class for use by other classes
#	in particular, creates 'tabLDA' for initialization
#	of the profile by LDA
#
###########################################################

class param(object):
    def __init__(self):
    	self.c = c
        self.V0 = V0
        self.V = V
        self.xspace = xspace
        self.dt = dt
        self.Delta_t = Delta_t
        self.tmax = tmax
        self.LDAtab = np.zeros((100,5)) #contains data needed for LDA
        								#each line is (thetaF,mu,rho,v,K)
        #for calculation of LDAtab
        varphi = lambda k1, k2: 2*self.c/((k1-k2)**2 + self.c**2)
        self.LDAtab[0,:] = np.array([0.,0.,0.,0.,1.])
        count=1
        for theta in np.linspace(0.001,5.,100-1):
        	lam_discr = np.linspace(-theta,theta,100)
        	varphi_mat = np.zeros((lam_discr.shape[0],lam_discr.shape[0]))
        	for j in range(lam_discr.shape[0]):
        		for i in range(lam_discr.shape[0]):
        			varphi_mat[i,j] = varphi( lam_discr[i], lam_discr[j] )
        	dlam = np.array([lam_discr[1]-lam_discr[0] for j in range(lam_discr.shape[0])])
        	dlam[0]=dlam[0]/2.; dlam[-1]=dlam[-1]/2.;
        	dresser = np.eye(lam_discr.shape[0]) - 1/(2*np.pi) * varphi_mat.dot(np.diag(dlam))
        	h0 = np.array( [1. for p in range(lam_discr.shape[0])] )
        	h1 = np.array( [lam_discr[p] for p in range(lam_discr.shape[0])] )
        	h2 = np.array( [lam_discr[p]**2/2. for p in range(lam_discr.shape[0])] )
        	h0dr = np.linalg.solve(dresser, h0)
        	h1dr = np.linalg.solve(dresser, h1)
        	h2dr = np.linalg.solve(dresser, h2)
        	self.LDAtab[count,0] = theta
        	self.LDAtab[count,1] = h2dr[0]/h0dr[0]
        	self.LDAtab[count,2] = 1/(2*np.pi) * h0dr.dot(dlam)
        	self.LDAtab[count,3] = h1dr[-1]/h0dr[-1]
        	self.LDAtab[count,4] = h1dr[0]**2.
        	count+=1
#        print self.LDAtab