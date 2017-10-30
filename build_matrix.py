# build_matrix.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# APPM 4380: Project 3. Least Squares Inversion
# Code to build equation matrix.
#
# The Algorithm works by counting the number of points along the trajectory
# of the X-Ray in each grid-box. The number of points tallied correspond
# to the coefficient in the equation matrix.
#
# For rays up to pi/4, the equation generates a list of points in x
# spanning -sqrt(2),sqrt(2), computes y = tan(theta) + C where C determines
# the initial vertical position of the ray. If the angle is greater than
# pi/4., the array generates a vector of y and uses x = y*cot(theta) + C_x
#

__author__ = "Omkar H. Ramachandran"

import numpy as np
import math

# Grid dimensions
N = 60
# Number of rays per angle
m = 10
# Number of points per ray
n = 120

# Generate Theta array
theta_i = 0.
theta_f = np.pi/2.
Ntheta = 100
th = np.linspace(theta_i,theta_f,Ntheta)

# Define shape of A matrix
A = np.zeros([m*Ntheta,N**2])

# Do the cot theta division and function call outside the loop.
# multiplication is ~ 40 times faster chip level than sin/cos/tan
tanth = np.tan(th)
cotth = 1./np.tan(th)

for t in range(Ntheta):
	x = np.linspace(-1,1,n)
	C = np.linspace(-1,1,m)
	y = np.zeros([m,n])
	for i in range(m):
		# Define straight line corresponding to the ray
		# Because y is a 2D array and x is 1D, to accomodate the floating
		# constant, I'm not changing names when switching from equation
		# in x for y (y=tan(theta)*x+C) and (x=cot(theta)*y+C). The reason
		# why this works is that when 0<th<pi/4, every point in x is 
		# sampled while for pi/4<th<pi/2, every point in y is sampled,
		# thus it is sufficient to define one spanning vector and define
		# the other in terms of it
		if(th[t] > np.pi/4.):
			y[i] = x*cotth[t]+C[i]
		else:
			y[i] = x*tanth[t]+C[i]

		for j in range(1,n-1):
			X = int((x[j]+1)/2.*(N))
			Y = int((y[i][j]+1)/2.*(N))
			print(X,Y,":",x[j],y[i,j],":",j,":",m)
			# Ignore points exactly on the boundary of the grid
			if(abs(x[j]) >= 1  or abs(y[i][j]) >= 1):
				continue
			# ... Otherwise, tally up that point in A
			A[t*m+i,N*X+Y] += 1
