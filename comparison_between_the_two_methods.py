#
# python code to compute points on a sphere
# uses a simple Monte Carlo method 
# draws the points on the sphere and plots the energy against iteration number

"""STOPPING CRITERIA: ONCE REACHES LITERATURE VALUE"""

import numpy as np
import matplotlib.pyplot as plt
import random
#from mpl_toolkits.mplot3d import Axes3D
import time

# project one or all points to the sphere
def initialproj(x):
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        while norm > 1.0:
            x[i] = (2.0*np.random.random(3)-1.0) #creates new random point if norm>1
            norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
    return x
    
def proj(x) :
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
    return x
    
def montecarlo(x,energy,loops,often,amplitude,ee):
    loop = 0
    while energy>ee:
        loop = loop + 1
        # randomly choose a point to move
        i=random.randint(0,n-1)
        # store the old coordinates of this point
        old=np.array(x[i])
        # randomly move this point
        x[i]=x[i]+amplitude*(2.0*np.random.random(3)-1.0)
        # put this point back on the sphere
        norm=np.sqrt(sum(x[i]**2))
        x[i]=x[i]/norm
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
        if(loop%often==0):
            looptime = time.time() -starttime
            #print("{0} {1:.6f}".format(loop,energy))
            out = open('out','a')
            out.write("{0} {1:.6f} \n".format(loop,energy))
            out.close()
            mclooplist.append(loop)
            mcenergylist.append(energy)
            mctimelist.append(looptime)
    return x, energy
    
def calcenergy(x):
    energyval=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            energyval=energyval+1.0/distance
    return energyval
    
def steepestdecent(x, amplitude): 
    forceslist = [] #empty list we will add the resultant components of each xi to
    forces_pointsjoni = np.zeros([n,3]) #empty array which we will add the resultant force of each point on a single other point
    for i in range(0,n):
        for j in range(0,n):
            if i==j:
                forces_pointsjoni[j] = 0
            else:
                forces_pointsjoni[j] = ((x[i]-x[j])/(np.sqrt(sum((x[i]-x[j])**2))**3)) #resultant forces of point j on point i
        resultanttotals = np.sum(forces_pointsjoni, axis=0) #sum up the columns to get total x,y,z compnents of forces acting on i
        forcetilda = resultanttotals - np.dot(resultanttotals, x[i]) * x[i] #using Gram-Schmidt to remove radial component
        forceslist.append(forcetilda) #adding our force-tilda to our empty list
    for i in range(0,n):        
        x[i] = x[i] + amplitude*forceslist[i] #moving all the points by the force-tilda multiplied by step size
    return x
    
def relax(x,amplitude,energy,ee):
    loop = 0
    while energy > ee:  
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
        
        looptime = time.time() -starttime
        if looptime == 0.00000000:
            looptime = 0.000001
            
        mosdlooplist.append(loop)
        mosdenergylist.append(energy)
        mosdtimelist.append(looptime)
        
        mosde = open('MOSDe.txt','a')
        mosde.write('%s %s \n' % (looptime, energy))
        mosde.close()
        
    return x, energy
 
                
                                              
# set the number of points
"""Choose either 11 or 12"""
#n=11
#ee=40.59646
n=12
ee = 49.16526





# set the number of loops
loops=10000
# set how often to output the energy 
often=100
# set the maximum amplitude of the change
amplitude=0.01

mctimes = []
mosdtimes = []

mosde = open('MOSDe.txt','w+')


#assign random start points on the sphere
random.seed()
initialx=initialproj((2.0*np.random.random((n,3))-1.0))

# calculate the initial energy
initialenergy=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance=np.sqrt(sum((initialx[i]-initialx[j])**2))
        initialenergy=initialenergy+1.0/distance
        
mosde = open('MOSDe.txt','a')
mosde.write('%s %s \n' % ('0.000001', initialenergy))
mosde.close()

x = initialx.copy()

"""MONTE CARLO"""
mclooplist=[0]
mcenergylist=[initialenergy]
mctimelist=[0.0000001]
starttime = time.time()
mcx, mcfinalenergy = montecarlo(initialx,initialenergy,loops,often,amplitude,ee)      
#print("Final energy, MC = {0:.6f} \n".format(mcfinalenergy)) 
endtime = time.time()
tt = endtime-starttime
mctimes.append(tt)

"""MOSD""" #need to rename different energies
mosdlooplist=[0]
mosdenergylist=[initialenergy]
mosdtimelist=[0.0000001]
starttime = time.time()
mosdx, mosdfinalenergy = relax(x,0.1,initialenergy,ee)
#print("Final energy, MOSD = {0:.6f} \n".format(mosdfinalenergy))   
endtime = time.time()
tt = endtime-starttime
mosdtimes.append(tt)


#print mosdenergylist[0]
#print mosdenergylist[1]

#print mosdtimelist[0]
#print mosdtimelist[1]

maxx1 = mosdtimelist[-1] + 0.2
maxx2 = mctimelist[-1] + 1.5

#Plot
fig1 = plt.figure(1)

ax4=fig1.add_subplot(212)
ax4.plot(mctimelist,mcenergylist, color='purple', label='Monte Carlo method')
ax4.plot(mosdtimelist,mosdenergylist, color='red', label='Method of steepest descent')
ax4.set_xlabel("time elapsed, seconds")
ax4.set_ylabel("energy")
ax4.set_title("Method of steepest descent and Monte Carlo method")
ax4.legend(loc='best', fontsize='11')
ax4.semilogx()

ax5=fig1.add_subplot(211)
ax5.plot(mctimelist,mcenergylist, color='purple', label='Monte Carlo method')
ax5.plot(mosdtimelist,mosdenergylist, color='red', label='Method of steepest descent')
ax5.set_xlabel("time elapsed, seconds")
ax5.set_ylabel("energy")
ax5.set_title("Method of steepest descent and Monte Carlo method")
ax5.legend(loc='upper right', fontsize='11')

plt.tight_layout(1)

plt.show() 



