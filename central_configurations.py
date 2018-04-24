# python code to compute points in a ball
# uses a method of steepest descent
# plots the points in a ball, with the points coloured according to their radius

import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
import time
 

starttime = time.time()


# project one or all points to the sphere
   
def calcenergy(x):
    energyval=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            energyval=energyval+1.0/distance
    for i in range(0,n):
        r_i = np.dot(x[i],x[i])
        energyval = energyval + 0.5*r_i
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
        resultanttotals = np.nansum(forces_pointsjoni, axis=0) - x[i] #sum up the columns of get total x,y,z compnents of forces acting on i
        forceslist.append(resultanttotals) #adding our forcetilda to our empty list
    for i in range(0,n):
        x[i] = x[i] + amplitude*forceslist[i]
    return x
    
def relax(x,amplitude,energy):
    loop = 0
    while amplitude > 0.00016:  
        loop = loop+1
        #store the value of the energy to compare to when we calucalte a new energy
        oldenergy = energy.copy()
        # store the old coordinates of this point
        oldx = np.array(x[i])
        #calculate resultant forces of each point
        x = steepestdecent(x, amplitude)
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
    return x, energy

n=24


   
k = 100    
# set the number of points
while k >4:
    
    looplist=[]
    energylist=[]
    
    #assign random start points on the sphere
    random.seed()
    x = (2*(n**(1/3))*np.random.random((n,3))-1*(n**(1/3)))
    
    # calculate the initial energy
    energy=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance = abs(sum(x[i]-x[j]))
            energy=energy+1.0/distance
    for i in range(0,n):
        r_i = np.dot(x[i],x[i])
        energy = energy + 0.5*r_i
        
    x, finalenergy = relax(x,0.5,energy)
    
    #find radii of different layers- 'distances' gives list of radii
    rilist = []
    for i in range(0,n):
        ri = np.sqrt(np.dot(x[i],x[i]))
        #print ri
        rilist.append(ri)
    distances = [ri]
    for i in range(0,len(rilist)):
        count = 0
        for j in range(0,len(distances)):
            diff = abs(rilist[i] - distances[j])
            if diff > 0.2:
                count = count + 1
                if count == len(distances):
                    if rilist[i] > 0.5:
                        distances.append(rilist[i])
    maxdist = 0
    for i in range(0,len(distances)):
        if distances[i] > maxdist:
            maxdist = distances[i]
    
    rings = len(distances)
    k = rings
    
    if k<4:
        # output final energy to the screen and points to a file
        print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy))
        
            

#Render
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')

for r in range(0,rings):

    #Create a sphere
    theta, phi = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
    xs = distances[r]*np.sin(theta)*np.cos(phi)
    ys = distances[r]*np.sin(theta)*np.sin(phi)
    zs = distances[r]*np.cos(theta)
    ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.3, linewidth=0)
    
#convert data
x1=[]
x2=[]
x3=[]
for i in range(0,n):
    x1.append(x[i,0])
    x2.append(x[i,1])
    x3.append(x[i,2])
    

for i in range(0,n):
    if abs(np.sqrt(x1[i]**2 + x2[i]**2 + x3[i]**2) - maxdist) < 0.2:
        ax.scatter(x1[i],x2[i],x3[i],s=80, c='blue')
    else:
        ax.scatter(x1[i],x2[i],x3[i],s=80, c='red')
        
plt.axis('off')
ax.view_init(elev=0.,azim=0)

ax.set_xlim([-maxdist,maxdist])
ax.set_ylim([-maxdist,maxdist])
ax.set_zlim([-maxdist,maxdist])
ax.set_aspect("equal")
ax.set_title("{0} ".format(n)+"points on a sphere")

plt.show() 

ppp = np.sqrt(sum((x[0])**2))
#print ppp

endtime = time.time()
#print(endtime-starttime)

