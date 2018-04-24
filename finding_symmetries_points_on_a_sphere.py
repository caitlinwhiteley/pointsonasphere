# python code to compute symmetries of points on a sphere
# uses a method of steepest descent
# plots the points on a sphere

import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D



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
                
# set the number of points
#n=int(raw_input('n?'))
n = 25

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
        

    # output energy to screen and a file
    if(loop%often==0):
        #print("{0} {1:.6f}".format(loop,energy))
        out = open('out','a')
        out.write("{0} {1:.6f} \n".format(loop,energy))
        out.close()
        looplist.append(loop)
        energylist.append(energy) 
        

# output final energy to the screen and points to a file
print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(energy))
        
points=open('points.out','w')
for i in range(0,n):
    for j in range(0,3):
        points.write("{0:.6f} ".format(x[i,j]))
    points.write('\n')              
points.close()


# find dipole moment
dipolec = np.sum(x, axis=0)
dipole = np.sqrt(dipolec[0]**2 + dipolec[1]**2 + dipolec[2]**2)

if dipole < 0.000005:
    dipole = 0



#working out moment of inertia tensor

ixx = 0
iyy = 0
izz = 0
ixy = 0
ixz = 0
iyz = 0

for i in range(0,n):
    nx = x[i,0] #column x1
    ny = x[i,1] #column x2
    nz = x[i,2] #column x3
    
    ny2nz2 = ny**2 + nz**2 #calculating x2 squared + x3 squared for row i (equation for inertia tensor)
    nx2nz2 = nx**2 + nz**2
    nx2ny2 = nx**2 + ny**2
    nxny = nx*ny
    nxnz = nx*nz
    nynz = ny*nz    
        
    ixx = ixx + ny2nz2 #adding the calculations for each row on
    iyy = iyy + nx2nz2
    izz = izz + nx2ny2
    ixy = ixy + nxny
    ixz = ixz + nxnz
    iyz = iyz + nynz
    
#putting all our calculated values into a matrix I
matrixi = np.array([[ixx, -ixy, -ixz], [-ixy, iyy, -iyz], [-ixz, -iyz, izz]])
    

#find eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(matrixi)

#print eigenvalues

sameeigens = 1

#checking for if there are 2 eigenvalues the same. If there are, we work out which eigenvalue is the
#one without a pair, as we will make the corresponding column in x swap with the 1st column in x
for i in range(0,2): 
    for j in range(i+1,3):
        absvalue = abs(eigenvalues[i] - eigenvalues[j])
        if absvalue < 0.001:
            sameeigens = 2
            if i == 0 and j == 1:
                zaxis = 2
            elif i == 0 and j == 2:
                zaxis = 1
            elif i == 1 and j == 2:
                zaxis = 0
            

#check if all 3 eigenvalues are the same                      
if abs(((eigenvalues[0]+eigenvalues[1]+eigenvalues[2])/3) - eigenvalues[0]) < 0.001:
    sameeigens = 3
    #print("platonic symmetry")
    

#copy our old x matrix before we start swapping columns around
oldx = x.copy()

#if we have 2 repeated eigenvalues, then we need to swap columns around so the single eigenvalue
#corresponds to the first column in x
if sameeigens == 2:
    if zaxis == 0: #if the single eigenvalue is already in the first column
        newX = x
    elif zaxis == 1:
        xcol = x[:,0].copy()
        ycol = x[:,1].copy()
        zcol = x[:,2].copy()
        intermediateX = np.array([ycol,xcol,zcol]) #this puts the 3 copied columns as 3 rows
        newX = intermediateX.transpose(1,0) #transpose to get 3 columns again
    elif zaxis == 2:
        xcol = x[:,0].copy()
        ycol = x[:,1].copy()
        zcol = x[:,2].copy()
        intermediateX = np.array([zcol,ycol,xcol])
        newX = intermediateX.transpose(1,0)
        
if sameeigens == 3 or sameeigens == 1:
    newX = x
    
  

rotateddata2 = np.zeros([n,3]) #set up empty array to put rotated data in

for i in range(0,n):
    rotateddata2[i] = np.dot(eigenvectors.T, newX[i]) #rotate x by multiplying with eigenvectors

xxcol = rotateddata2[:,0].copy()
yycol = rotateddata2[:,1].copy()
zzcol = rotateddata2[:,2].copy()
rotateddata3 = np.array([yycol,zzcol,xxcol])
rotateddata = rotateddata3.T
#rotateddata = rotateddata2
 

def compare(xcol,ycol,zcol,newx,newy,newz):
    tallie = 0
    for i in range(0,n):
        mindistance = 5
        for j in range(0,n):
            xdist = xcol[i] - newx[j]
            ydist = ycol[i] - newy[j]
            zdist = zcol[i] - newz[j]
            distance = np.sqrt(xdist**2 + ydist**2 + zdist**2)
            if distance < mindistance:
                mindistance = distance
        if mindistance < 0.1:
            tallie = tallie + 1
        else:
            break
    return tallie
            
cnsymmetry = 100
coords = rotateddata.T
xcol, ycol, zcol = np.hsplit(rotateddata,3)

if sameeigens == 2: #NEEDS TO BE C3 SYMMETRY
    for mm in range(0,12):
        m = 12 - mm #so we count down from 12
        alpha = 2*np.pi/m
        rotmatrix = np.array([[np.cos(alpha),-np.sin(alpha),0],[np.sin(alpha),np.cos(alpha),0],[0,0,1]])        
        symmetrym = np.dot(rotmatrix, coords)
        asd = symmetrym.T
        newx, newy, newz = np.hsplit(asd,3)
        
        if compare(xcol,ycol,zcol,newx,newy,newz) == n:
            cnsymmetry = m
            break


       
if sameeigens == 1:
    symm = "at most"
    k = "C 2"
elif dipole > 0:
    symm = "C"
    k = cnsymmetry
elif sameeigens == 3:
    symm = "platonic"
    k = ""
elif dipole == 0:
    symm = "D"
    k = cnsymmetry

print "N =", n, ", ", symm, k, "symmetry"

if dipole >0:
    print "Dipole = ", dipole       



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
    x1.append(rotateddata2[i,0])
    x2.append(rotateddata2[i,1])
    x3.append(rotateddata2[i,2])
    

#Render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xs, ys, zs,  rstride=4, cstride=4, color='yellow', alpha=0.3, linewidth=0)
ax.scatter(x1,x2,x3,color="black",s=80)
plt.axis('off')
ax.view_init(elev=0.,azim=0)

ax.set_xlim([-1.0,1.0])
ax.set_ylim([-1.0,1.0])
ax.set_zlim([-1.0,1.0])
ax.set_aspect("equal")
ax.set_title("{0} ".format(n)+"points on a sphere")

plt.show() 



