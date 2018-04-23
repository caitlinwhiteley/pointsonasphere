# python code to compute points on a sphere
# uses a method of steepest descent
# draws the points on the sphere and plots the energy against iteration number


import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import time
 

starttime = time.time()


# project one or all points to the sphere
def proj(x) :
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
    return x
    
def initialproj(x):
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        while norm > 1.0:
            x[i] = (2.0*np.random.random(3)-1.0) #creates new random point if norm>1
            norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
    return x
    
def calcenergy(x):
    energyval=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            energyval=energyval+1.0/distance
    return energyval
    
def steepestdecent(x, amplitude): #called within relax function
    forceslist = [] #empty list we will add the resultant components of each xi to
    forces_pointsjoni = np.zeros([n,3]) #empty array which we will add the resultant force of each point on a single other point
    for i in range(0,n):
        for j in range(0,n):
            if x[i,2:]==x[j,2:]:
                forces_pointsjoni[j] = 0
            else:
                forces_pointsjoni[j] = ((x[i]-x[j])/(np.sqrt(sum((x[i]-x[j])**2))**3)) #resultant forces of point j on point i
        resultanttotals = np.nansum(forces_pointsjoni, axis=0) #sum up the columns of get total x,y,z compnents of forces acting on i
        forcetilda = resultanttotals - np.dot(resultanttotals, x[i]) * x[i] #taking dot product to only take forces in relevant direction
        forceslist.append(forcetilda) #adding our forcetilda to our empty list
    for i in range(0,n):
        x[i] = x[i] + amplitude*forceslist[i]
    return x
    
def relax(x,amplitude,energy):
    loop = 0
    while amplitude > 0.0016:  
        loop = loop+1
        #store the value of the energy to compare to when we calucalte a new energy
        oldenergy = energy.copy()
        # store the old coordinates of this point
        oldx = np.array(x[i])
        #calculate resultant forces of each point
        x = steepestdecent(x, amplitude)
        x = proj(x)
        # calculate the energy
        energy = calcenergy(x)        
        #reducing step size if energy increases 
        difference = energy-oldenergy
        if(difference>0.0): #if the old energy is lower than the new energy reduce step size
            amplitude = amplitude/2
            x[i] = oldx  
            energy = oldenergy
        #calcualte energy after iterations 
        looplist.append(loop)
        energylist.append(energy)   
    return x, energy, loop
    
                
# set the number of points
n = 11

looplist=[]
energylist=[]

#assign random start points on the sphere
random.seed()
x=initialproj((2.0*np.random.random((n,3))-1.0))

# calculate the initial energy
energy=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance = abs(sum(x[i]-x[j]))
        energy=energy+1.0/distance

x, finalenergy, loops = relax(x,0.2,energy)


# output final energy to the screen and points to a file
print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy))
  
#print loops


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
ax = fig.add_subplot(211, projection='3d')
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.3, linewidth=0)
ax.scatter(x1,x2,x3,s=80, color='black')
plt.axis('off')
ax.view_init(elev=0.,azim=0)

ax.set_xlim([-1.0,1.0])
ax.set_ylim([-1.0,1.0])
ax.set_zlim([-1.0,1.0])
ax.set_aspect("equal")
ax.set_title("{0} ".format(n)+"points on a sphere")

ax2=fig.add_subplot(212)
ax2.plot(looplist,energylist)
ax2.set_xlabel("iteration number")
ax2.set_ylabel("energy")

plt.show() 

endtime = time.time()
#print(endtime-starttime)

