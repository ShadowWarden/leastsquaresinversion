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
import scipy.sparse as scp
import math

RHO0 = 3

def update(A,X,Y,N,th,b,S):
    if(X[i] >= N or Y[i] >= N or Y[i] < 0 or X[i] < 0):
        return
    if(abs(th) < np.pi/4):
        A[N*Y[i]+X[i]] += 1
        b += S[2*Y[i],2*X[i]]
    else:
        A[N*X[i]+Y[i]] += 1
        b += S[2*X[i],2*Y[i]]
    return A, b
# Grid dimensions
N = 20
# Number of rays per angle
m = 2*N
# Number of points per ray
n = 2*N

# Generate Theta array
theta_i = -np.pi/2.
theta_f = np.pi/2.
Ntheta = 120
th = np.linspace(theta_i,theta_f,Ntheta)

# Define shape of A matrix
A = np.zeros((m*Ntheta,N**2))
b = np.zeros([m*Ntheta])

xx = np.linspace(-1,1,N)
yy = np.linspace(-1,1,N)

Xx,Yy = np.meshgrid(xx,yy)

S = RHO0*np.heaviside(1.-Xx**2-Yy**2,0.5)*(0.7-Xx**2-Yy**2)

#S = np.genfromtxt("CU_logo.csv",delimiter=",")

# Do the cot theta division and function call outside the loop.
# multiplication is ~ 40 times faster chip level than sin/cos/tan
        
for t in range(Ntheta):
        x = np.linspace(-np.sqrt(2),np.sqrt(2),n)
        C = np.linspace(-np.sqrt(2),np.sqrt(2),m)
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
            if(abs(th[t]) > np.pi/4.):
                y[i] = -abs(th[t])/th[t]*(x*np.tan(np.pi/2 - abs(th[t])) -C[i])
            else:
                y[i] = -abs(th[t])/th[t]*(x*np.tan(th[t]) - C[i])    
            for j in range(n):
                X = (x[j]+1)/2.*(N)
                Y = (y[i][j]+1)/2.*(N)
            # Ignore points exactly on the boundary of the grid
            # ii = np.unique(np.append(ii,jj))
            # X = np.delete(X,ii)
            # Y = np.delete(Y,ii)
            # X[ii] = N-1
            # Y[jj] = N-1
                X = X.astype(int)
                Y = Y.astype(int)
            # ... Otherwise, tally up that point in A. Vectorize for
            # maximum performance
                if(X < 0 or Y < 0 or X >= N or Y >= N):
                    continue
                if(abs(th[t]) < np.pi/4):
                    A[t*m+i,N*Y+X] += 1
                    b[t*m+i] += S[Y,X]
                else:
                    A[t*m+i,N*X+Y] += 1
                    b[t*m+i] += S[X,Y]

ii = np.where(np.all(A == 0, axis=1))
A = A[~np.all(A == 0, axis=1)]
b = np.delete(b,ii)
#for i in range(len(A)):
#    A[i] /= np.linalg.linalg.norm(A[i])
#    b[i] /= np.linalg.linalg.norm(A[i])
