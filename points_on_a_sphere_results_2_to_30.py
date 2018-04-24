# python code to compute points on a sphere using a genetic algorithm
# uses a method of steepest descent

import numpy as np
#import matplotlib.pyplot as plt
import random
#'from mpl_toolkits.mplot3d import Axes3D


# project one or all points to the sphere
def proj(x) :
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
    return x
    
def initialproj(x):
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        if norm > 1.0:
            x[i] = (2.0*np.random.random(3)-1.0) #creates new random point if norm>1
            norm=np.sqrt(sum(x[i]**2))
            x[i]=x[i]/norm
        else:
            x[i]=x[i]/norm
    return x

points = open('Points.txt','w+')
energies = open('Energies.txt','w+')

listn = []
listenergies = []

for q in range(2,31):
                    
    # set the number of points
    #n=int(raw_input('n?'))
    n = q
    
    # open file to output energy during minimization 
    out = open('out','w')
    out.close()
    
    looplist=[]
    energylist=[]
    
    #assign random start points on the sphere
    random.seed()
    x=initialproj((2.0*np.random.random((n,3))-1.0))
    
    
    # calculate the initial energy
    energy=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            energy=energy+1.0/distance
    
    
    amplitude = 0.2
    
    # the main loop to reduce the energy        
    loop = 0
    
    while amplitude > 0.0016:
        
        loop = loop + 1
            
        #store the value of the energy to compare to when we calucalte a new energy
        oldenergy = energy.copy()
    
        # store the old coordinates of this point
        oldx = np.array(x[i])
    
    
        #calculate resultant forces of each point
        
        forceslist = [] #empty list we will add the resultant components of each xi to
        forces_pointsjoni = np.zeros([n,3]) #empty array which we will add the resultant force of each point on a single other point
    
        for i in range(0,n):
            for j in range(0,n):
                if i == j:
                    forces_pointsjoni[j] = 0
                else: forces_pointsjoni[j] = ((x[i]-x[j])/(np.sqrt(sum((x[i]-x[j])**2))**3)) #resultant forces of point j on point i
            resultanttotals = np.nansum(forces_pointsjoni, axis=0) #sum up the columns of get total x,y,z compnents of forces acting on i
            forcetilda = resultanttotals - np.dot(resultanttotals, x[i]) * x[i] #taking dot product to only take forces in relevant direction
            forceslist.append(forcetilda) #adding our forcetilda to our empty list
    
            x[i] = x[i] + amplitude*forceslist[i]
            
        x = proj(x)
            
        
        # calculate the energy
        energy=0.0
        for i in range(0,n):
            for j in range(i+1,n):
                distance=np.sqrt(sum((x[i]-x[j])**2))
                energy=energy+1.0/distance
                
    
        #reducing step size if energy increases 
        difference = energy-oldenergy
        if(difference>0.0): #if the old energy is lower than the new energy reduce step size
            amplitude = amplitude/2
            x[i] = oldx
            
    
    # output final energy to the screen and points to a file
    print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(energy))
        

    listn.append(q)
    listenergies.append(energy)   
    
    points.write('%s \n' % n)
    for i in range(0,n): points.write('%s %s %s \n' % (x[i,0],x[i,1],x[i,2]))
    energies.write('%s %s \n' % (n, energy))

 
points.close()
energies.close()    
