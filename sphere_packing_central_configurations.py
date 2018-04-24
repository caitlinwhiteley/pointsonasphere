# python code to compute points in a ball, using sphere packing representation
# uses a method of steepest descent


import numpy as np
import matplotlib.pyplot as plt
import random
#from mpl_toolkits.mplot3d import Axes3D
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
            if i==j:
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
    
                
# set the number of points
n = 11

looplist=[]
energylist=[]

#assign random start points on the sphere
random.seed()
x = (2.0*(n**(1/3))*np.random.random((n,3))-(n**(1/3)))

# calculate the initial energy
energy=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance = abs(sum(x[i]-x[j]))
        energy=energy+1.0/distance
for i in range(0,n):
    r_i = np.dot(x[i],x[i])
    energy = energy + 0.5*r_i
    
x, finalenergy = relax(x,0.3,energy)






# output final energy to the screen and points to a file
print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy))
  

#find radii of different layers- 'distances' gives list of radii
rilist = []
for i in range(0,n): #list of all the radii
    ri = np.sqrt(np.dot(x[i],x[i])) #calculate distance of point i from origin
    rilist.append(ri)
distances = [ri]
for i in range(0,len(rilist)):
    count = 0
    for j in range(0,len(distances)):
        diff = abs(rilist[i] - distances[j])
        if diff > 0.1:
            count = count + 1
            if count == len(distances):
                distances.append(rilist[i])

maxdist = 0
for i in range(0,len(distances)):
    if distances[i] > maxdist:
        maxdist = distances[i]

rings = len(distances)

#print rilist

#find min distance between points = R
mindistt = 10000
for i in range(0,n):
    for j in range(i+1,n):
        distt = np.sqrt(sum((x[i]-x[j])**2))
        if distt < mindistt:
            mindistt = distt
minr = mindistt/2

#Render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#sphere packing all same colour. For different colours, make lists            
for i in range(0,n): #range(0,len(a))
    u = np.linspace(0,  2*np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    xs = minr * np.outer(np.cos(u), np.sin(v)) + x[i,0]
    ys = minr * np.outer(np.sin(u), np.sin(v)) + x[i,1]
    zs = minr * np.outer(np.ones(np.size(u)), np.cos(v)) + x[i,2]
    if len(distances) == 1:
        ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='plum', alpha=0.8, linewidth=0)
    if len(distances) == 2: #here alpha determines hue, not transparency
        if abs(np.sqrt(np.dot(x[i],x[i]))-distances[0]) < 0.1:
            ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='blue', alpha=0.8, linewidth=0)
        else:
            ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='plum', alpha=0.8, linewidth=0)
    if len(distances) == 2:
        if abs(np.sqrt(np.dot(x[i],x[i]))-distances[0]) < 0.1:
            ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='green', alpha=0.9, linewidth=0)
        if abs(np.sqrt(np.dot(x[i],x[i]))-distances[1]) < 0.1:
            ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='blue', alpha=0.8, linewidth=0)
        else:
            ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='plum', alpha=0.8, linewidth=0)
     

plt.axis('off')
ax.view_init(elev=0.,azim=0)

ax.set_xlim([-maxdist,maxdist])
ax.set_ylim([-maxdist,maxdist])
ax.set_zlim([-maxdist,maxdist])
ax.set_aspect("equal")
ax.set_title("{0} ".format(n)+"points on a sphere")

plt.show() 

endtime = time.time()
#print(endtime-starttime)

