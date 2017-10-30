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
m = 100
# Number of points per ray
n = 20

theta_i = 0
theta_f = np.pi/4.
Ntheta = 1
th = np.linspace(theta_i,theta_f,Ntheta)

A = np.zeros([m*Ntheta,N**2])

for t in range(Ntheta):
	x = np.linspace(-np.sqrt(2),np.sqrt(2),n)
	C = np.linspace(-np.sqrt(2),np.sqrt(2),m)
	y = np.zeros([m,n])
	for i in range(m):
		y[i] = x*np.tan(th[t])+C[i]
		for j in range(n):
			if((x[j] >= -1 and x[j] <= 1) and (y[i][j] >= -1 and y[i][j] <= 1)):
				X = int((x[j]+1)/2./np.sqrt(2)*(N))
				Y = int((y[i][j]+1)/2./np.sqrt(2)*(N))
				print(X,Y,":",x[j],y[i,j])
				A[t*Ntheta+i,N*Y+X] += 1
