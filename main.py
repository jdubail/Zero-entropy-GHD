""" Calculate the zero-entropy GHD solution, for parameters """
""" defined in file 'parameters.py'. 						"""
""" The files 'cont_tX.dat' contain the contour in			"""
""" phase-space at time X, and the files 'density_tX.dat'	"""
""" contain the particle density.							"""
""" Based on PRL 119, 195301 (2017)							"""
""" (https://doi.org/10.1103/PhysRevLett.119.195301) 		"""


#libraries
import numpy as np
import os

import parameters
import TBA
import contour

#initialize parameters according to file 'parameters.py'
para = parameters.param()
#initialize contour
contour = contour.contour(para)


#create folder 'data' if it does not exist
if not os.path.exists('data'):
    os.mkdir('data')

#main loop for time evolution
for p in range(int(para.tmax/para.dt)+2):
	if p%int(para.Delta_t/para.dt)==0:	#save files
		print 'save'
		np.savetxt('data/cont_t%d.dat' % (p/int(para.Delta_t/para.dt)), contour.pts)
		np.savetxt('data/density_t%d.dat' % (p/int(para.Delta_t/para.dt)), contour.density())
	print 't =', p*para.dt
	contour.evolve(para.dt,para.V)