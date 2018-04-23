# python code to compute points on a sphere
# uses a simple Monte Carlo method 
# draws the points on the sphere 

import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import time

starttime = time.time()

# project one or all points to the sphere
def proj(x,j):
    if(j==n):
        for i in range(0,n):
            norm=np.sqrt(sum(x[i]**2))
            x[i]=x[i]/norm
    else:
        norm=np.sqrt(sum(x[j]**2))
        x[j]=x[j]/norm
    return x
    
def montecarlo(x,energy,loops,often,amplitude):
    for loop in range(0,loops+1):
        # randomly choose a point to move
        i=random.randint(0,n-1)
        # store the old coordinates of this point
        old=np.array(x[i])
        # randomly move this point
        x[i]=x[i]+amplitude*(2.0*np.random.random(3)-1.0)
        x=proj(x,i)
        # calculate the difference in energy
        difference=0.0
        for j in range(0,n):
            if(j!=i):
                distance=np.sqrt(sum((x[i]-x[j])**2))
                distanceold=np.sqrt(sum((old-x[j])**2))
                difference=difference+1.0/distance-1.0/distanceold;
        # accept or reject the move 
        if(difference<0.0):
            energy=energy+difference
        else:
            x[i]=old
    return x, energy
    
looplist=[]
energylist=[]
                
# set the number of points
n=12
# set the number of loops
loops=10000
# set how often to output the energy 
often=100
# set the maximum amplitude of the change
amplitude=0.01
#assign random start points on the sphere
random.seed()
x=proj((2.0*np.random.random((n,3))-1.0),n)

 
x=proj(x,n)

# calculate the initial energy
energy=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance=np.sqrt(sum((x[i]-x[j])**2))
        energy=energy+1.0/distance
     
x, finalenergy = montecarlo(x,energy,loops,often,amplitude)      
  
print("Final energy = {0:.6f} \n".format(finalenergy))      
        
    

#Create a sphere
theta, phi = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
xs = np.sin(theta)*np.cos(phi)
ys = np.sin(theta)*np.sin(phi)
zs = np.cos(theta)

#convert data
x1=[]
x2=[]
x3=[]
for i in range(0,n):
    x1.append(x[i,0])
    x2.append(x[i,1])
    x3.append(x[i,2])
    
#Render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.3, linewidth=0)
ax.scatter(x1,x2,x3,s=80, color='black')
plt.axis('off')
ax.view_init(elev=0.,azim=0)

ax.set_xlim([-1.0,1.0])
ax.set_ylim([-1.0,1.0])
ax.set_zlim([-1.0,1.0])
ax.set_aspect("equal")
ax.set_title("{0} ".format(n)+"points on a sphere")

plt.show() 


endtime = time.time()
#print(endtime-starttime)
