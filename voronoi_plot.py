# python code to compute points on a sphere and draw Voronoi plot
# uses a method of steepest descent 

import numpy as np
import matplotlib.pyplot as plt
import random


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
                

n = 14


# set how often to output the energy 
often=100

# set the maximum amplitude of the change
amplitude=0.05
 
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



       
s = 300

phi, theta = np.linspace(0,2*np.pi,s),np.linspace(0,np.pi,s)  
alpha,beta = np.meshgrid(phi,theta) 
xs = np.sin(beta)*np.cos(alpha) 
ys = np.sin(beta)*np.sin(alpha) 
zs = np.cos(beta) 
       
xxa = xs.flatten() 
yya = ys.flatten() 
zza = zs.flatten() 

nop=len(xxa)

xx = np.zeros([nop,4])

for i in range(0,nop):
    mindistance = 10
    for j in range(0,n):
        distancex = xxa[i] - x[j,0]
        distancey = yya[i] - x[j,1]
        distancez = zza[i] - x[j,2]
        distance = np.sqrt(distancex**2 + distancey**2 + distancez**2)
        if distance < mindistance:
            mindistance = distance
            colour = j
            xx[i] = np.array([xxa[i],yya[i],zza[i],j])
           

#Create a sphere
theta, phi = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
xs = np.sin(theta)*np.cos(phi)
ys = np.sin(theta)*np.sin(phi)
zs = np.cos(theta)

#convert data
xx1=[]
xx2=[]
xx3=[]
clist=[]
for i in range(0,nop):
    xx1.append(xx[i,0])
    xx2.append(xx[i,1])
    xx3.append(xx[i,2]) 
    clist.append(xx[i,3])    
         
                              
#Render
fig = plt.figure()
ax2 = fig.add_subplot(111, projection='3d')  
ax2.scatter(xx1,xx2,xx3,c=clist,edgecolor='none',cmap=plt.cm.plasma, s=5)

plt.axis('off')
ax2.view_init(elev=90.,azim=0)
            
ax2.set_xlim([-1.0,1.0])
ax2.set_ylim([-1.0,1.0])
ax2.set_zlim([-1.0,1.0])
ax2.set_aspect("equal")
ax2.set_title("points on a sphere")

plt.show() 

 



