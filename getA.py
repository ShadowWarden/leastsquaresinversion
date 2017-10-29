import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos,pi,sqrt
import pandas as pd

class A:
	def __init__(self,nodes,angles,rays):
		self.n = nodes
		self.m = angles
		self.nr = rays
		self.A = np.zeros([n*m,n*n])
		
		#make grid
		self.bounds = np.linspace(-n/2,n/2,n+1)
		
		self.grd_x,self.grd_y = np.meshgrid(self.bounds,self.bounds)

		
		length = len(self.bounds)
		self.arr_x =[]
		self.arr_y = []
		y_0 = np.array(range(n))
		x_0 = np.arange(0,n*n,n)
		for i in range(length-1):
			y_vec = y_0+i*n
			x_vec = x_0 + i
			self.arr_x.append(list(x_vec))
			self.arr_y.append(list(y_vec))
		
	
	def discretize_rays(self,theta):
		rotation_M = [[sin(theta),cos(theta)],[-sin(theta),cos(theta)]]
		X = np.asarray(self.grd_x).ravel()
		Y = np.asarray(self.grd_y).ravel()
		M = [list(X),list(Y)]
		MR = np.matmul(rotation_M,M)
		
		lower = min(MR[1])
		upper = max(MR[1])
		
		rays = np.linspace(lower,upper,self.nr)#[1:self.nr+1]
		
		discretized = np.linspace(sqrt(2)*min(MR[0]),sqrt(2)*max(MR[0]),self.n**3)#[1:self.n**3-1]
		return(rays,discretized)
		
		
	def pointToBox(self,x,y,theta,verbose=False):
		th = 2*pi-theta
		
		#invert rotation by rotating by 2pi-theta
		irot_M = [[cos(th),sin(th)],[-sin(th),cos(th)]]

		X,Y = np.matmul(irot_M,[x,y])

		X_ind = np.where(X <= self.bounds)[0]
		Y_ind = np.where(Y <= self.bounds)[0]
		
		if verbose:
			print(X_ind)
			print(Y_ind)
			
		if len(X_ind)==0 or len(Y_ind)==0:
			box = None
		else:
			index_X = min(X_ind)
			index_Y = min(Y_ind)
			if (index_X == 0 and not X==min(self.bounds)) or (index_Y==0 and not Y == min(self.bounds)):  #out of bounds of grid 
				box = None
			else:
				X_list = self.arr_x[index_X-1]
				Y_list = self.arr_y[index_Y-1]

				box = list(set(X_list)&set(Y_list))[0]
		
		return(box)
		
	def scan_data(self):
		
		angles = np.linspace(0,pi,self.m+1)[0:self.m]
		i=0
		rays,discretized = self.discretize_rays(0)
		
		j=1
		
		for theta in angles:
			print("Run "+str(j)+ " of "+str(self.m))
			j+=1
			for ray in rays:
				boxes = [0]*(n*n)
				for x in discretized:
					box = self.pointToBox(x,ray,theta)
					try:
						boxes[box]+=1
					except(TypeError):
						#print([x,ray])
						#self.pointToBox(x,ray,theta,verbose=True)
						pass
				maximum = max(boxes)
				self.A[i,:] = [a/maximum for a in boxes]
				i+=1
		print(self.A)
		
		df = pd.DataFrame(self.A)
		df.to_csv("A.csv",header=False,index=False)


if __name__ == "__main__":
	n = 62
	m = 66
	nr = 60
	
	scan = A(n,m,nr)
	scan.scan_data()

	