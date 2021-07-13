
""" class for zero-entropy GHD evolution of contour,	"""
""" based PRL 119, 195301 (2017)  (https://doi.org/10.1103/PhysRevLett.119.195301) """

#libraries
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d, splprep, splev

import TBA

#-------------------------------------------------------------------------
#	parameters
#-------------------------------------------------------------------------

N = 100		#number of points along the contour
ratio = 2.  #tunable parameter for uniformization of points along
			#the contour at time t=0

hbar = 1.05457e-25	#Planck's constant		(um^2.kg.ms^-1)
m = 1.44e-25		#mass    				(kg)

#-------------------------------------------------------------------------
#	class contour
#-------------------------------------------------------------------------

class contour(object):
    def __init__(self, param, filepts=None):	#initialize contour according to ground state
    											#of potential 'V0' in 'param.py'
		self.N = N								#number of points along the contour
		if filepts is None:
			self.pts = np.zeros((self.N,2))
			L=param.xspace.shape[0]
			fermi_x = np.zeros((param.xspace.shape[0]))
			v_x = np.zeros((param.xspace.shape[0]))
			theta_mu = interp1d(param.LDAtab[:,1],param.LDAtab[:,0], kind='linear')
			v_mu = interp1d(param.LDAtab[:,1],param.LDAtab[:,3], kind='linear')
			u1=0; u2=0;
			for u in range(L):
				if param.V0(param.xspace[u])<0:
					fermi_x[u] = theta_mu(-param.V0(param.xspace[u]))
					v_x[u] = v_mu(-param.V0(param.xspace[u]))
				if u>1 and param.V0(param.xspace[u])<0 and param.V0(param.xspace[u-1])>0: u1=u
				if u<L-1 and param.V0(param.xspace[u])<0 and param.V0(param.xspace[u+1])>0: u2=u
			curvey = np.zeros(2*(u2+1-u1)+1)
			curvex = np.zeros(2*(u2+1-u1)+1)
			for u in range(curvex.shape[0]/2):
				curvey[u] = (fermi_x[u1:u2+1])[u]
				curvex[u] = (param.xspace[u1:u2+1])[u]
			for u in range(curvex.shape[0]/2,curvex.shape[0]-1):
				curvey[u] = -(fermi_x[u1:u2+1])[curvex.shape[0]/2-1-u]
				curvex[u] = (param.xspace[u1:u2+1])[curvex.shape[0]/2-1-u]
			curvex[-1]=curvex[0]; curvey[-1]=curvey[0];
			tck, u = splprep([curvex/(max(curvex)-min(curvex)), curvey/(max(curvey)-min(curvey))/ratio], s=0)
			new_points = splev(np.linspace(0, 1, N+1), tck)
			new_points[0] = new_points[0][:N]*(max(curvex)-min(curvex))
			new_points[1] = new_points[1][:N]*(max(curvey)-min(curvey))*ratio
			self.pts[:,0]=new_points[0][:N]
			self.pts[:,1]=new_points[1][:N]
		else:
			ptsf=np.loadtxt(filepts)
			self.pts = ptsf
			if ptsf.shape[0]!=self.N:
				print "Warning, number of contour points is %d, not %d." % (ptsf.shape[0],self.N)
				self.N = ptsf.shape[0]	
    def evolve(self,dt,V):	#evolve one time step dt in potential V
		velocities = np.zeros(self.pts.shape[0])
		epsilon = 0.0000001	#for calculation of derivative
		#calculate effective velocities
		for j in range(self.pts.shape[0]):
			st = find_st(self.pts,j)
			ind = st.index(self.pts[j,1])
			sfs = TBA.SFS(st)
			sfs.calculate_veff()
			velocities[j] = sfs.veff[ind]
		#update pts on contour
		for j in range(self.pts.shape[0]):
			self.pts[j,1] += -(V(self.pts[j,0]+epsilon)-V(self.pts[j,0]))/epsilon * dt * hbar/m
			self.pts[j,0] += velocities[j] * dt * hbar/m
    def density(self):	#returns array with positions and densities
		dens = []
		for j in range(self.pts.shape[0]):
			st = find_st(self.pts,j)
			sfs = TBA.SFS(st)
			dens.append( [self.pts[j,0], sfs.charge_density(lambda u: 1.)] )
		dens.sort(key=lambda a: a[0])
		return np.array(dens)

#-------------------------------------------------------------------------
#	basic functions to deal with contour
#-------------------------------------------------------------------------

def find_pts(pts, j):	#given a point j, find the other segments [p-1,p]
						#at same position x 
	listp = []
	for p in range(pts.shape[0]):
		if pts[p-1,0]<pts[p,0]:
			if (pts[p-1,0]<=pts[j,0] and pts[j,0]<pts[p,0]): listp.append(p-1)
		else:
			if (pts[p,0]<pts[j,0] and pts[j,0]<=pts[p-1,0]): listp.append(p-1)
	return listp

def turningpt(pts,j):	#return 'True' if j is a turning point, otherwise 'False'
	if (pts[j,0]-pts[j-1,0])*(pts[j,0]-pts[(j+1)%pts.shape[0],0])>0: return True
	else: return False

def find_st(pts,j):		#finds the set of Fermi rapidities at point j
	listp = find_pts(pts, j)
	st = []
	for p in listp:
		if turningpt(pts,p)==False or ((p-j)%pts.shape[0])!=0:
			st.append( (pts[j,0]-pts[p,0])/(pts[p+1,0]-pts[p,0])*pts[p+1,1] \
				+(pts[p+1,0]-pts[j,0])/(pts[p+1,0]-pts[p,0])*pts[p,1] )
		else:
			st.append(pts[p,1])
			st.append(pts[p,1])
	st.sort()
	return st
