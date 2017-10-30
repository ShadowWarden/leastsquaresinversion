# build_matrix.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# APPM 4380: Project 3. Least Squares Inversion
# Code to build equation matrix
#

import numpy as np
import math

# Grid dimensions
N = 60
# Number of rays per angle
m = 10
# Number of points per ray
n = 20

# Generate Theta array
theta_i = 0.
theta_f = np.pi/2.
Ntheta = 100
th = np.linspace(theta_i,theta_f,Ntheta)

# Define shape of A matrix
A = np.zeros([m*Ntheta,N**2])

for t in range(Ntheta):
	# Define straight line corresponding to the ray
	x = np.linspace(-np.sqrt(2),np.sqrt(2),n)
	C = np.linspace(-np.sqrt(2),np.sqrt(2),m)
	y = np.zeros([m,n])
	for i in range(m):
		if(th[t] > np.pi/4.):
			y[i] = x/np.tan(th[t])+C[i]
		else:
			y[i] = x*np.tan(th[t])+C[i]

		for j in range(n):
			if((x[j] >= -1 and x[j] <= 1) and (y[i][j] >= -1 and y[i][j] <= 1)):
				X = int((x[j]+1)/2.*(N))
				Y = int((y[i][j]+1)/2.*(N))
				print(X,Y,":",x[j],y[i,j],":",j,":",m)
				if(th[t] > np.pi/4.):
					A[t*m+i,N*X+Y] += 1
				else:
					A[t*m+i,N*X+Y] += 1
