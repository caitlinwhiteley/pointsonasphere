# python code to compute points in a disc using a genetic algorithm
# uses a method of steepest descent
# does this for N=17,30,38,56,69
# these points are saved to a text file, which I then use in my code for points in a disc


import numpy as np
#import matplotlib.pyplot as plt
import random
#from mpl_toolkits.mplot3d import Axes3D
import time
#import sys 

starttime = time.time()


# project one or all points to the sphere
def proj(x) :
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        if norm>1:
            x[i]=x[i]/norm
    return x
    
def initialproj(x):
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        while norm > 1.0:
            x[i] = (2.0*np.random.random(2)-1.0) #creates new random point if norm>1
            norm=np.sqrt(sum(x[i]**2))
    return x
    
def calcenergy(x):
    energyval=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            if distance > 0:
                energyval=energyval+1.0/distance
    return energyval
    
def steepestdecent(x, amplitude): #called within relax function
    forceslist = [] #empty list we will add the resultant components of each xi to
    forces_pointsjoni = np.zeros([n,2]) #empty array which we will add the resultant force of each point on a single other point
    for i in range(0,n):
        for j in range(0,n):
            if i==j:
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
    while amplitude > 0.0003: 
        loop = loop+1 
        #store the value of the energy to compare to when we calucalte a new energy
        oldenergy = energy.copy()
        # store the old coordinates of this point
        oldx = x.copy()
        #calculate resultant forces of each point
        x = steepestdecent(x, amplitude)
        x = proj(x)
        # calculate the energy
        energy = calcenergy(x)        
        #reducing step size if energy increases 
        difference = energy-oldenergy
        if(difference>0.0): #if the old energy is lower than the new energy reduce step size
            amplitude = amplitude/2.0
            x = oldx.copy()
            #print difference
            #print loop
            energy = oldenergy.copy()
        #calcualte energy after iterations 
        looplist.append(loop)
        energylist.append(energy)      
    return x, energy
    
def bubblesort(newlist,xlist):
    for ii in range(len(newlist)-1,0,-1):
        for i in range(ii):
            if newlist[i]>newlist[i+1]:
                temp = newlist[i]
                newlist[i] = newlist[i+1]
                newlist[i+1] = temp
                temp3 = xlist[i]
                xlist[i] = xlist[i+1]
                xlist[i+1] = temp3
    return newlist, xlist
    
def randomrotation(x):
    random.seed()
    #choose random alpha, beta
    alpha = 2.0*np.pi*random.random()
    beta = 2.0*np.pi*random.random()
    #create 2 matrices for rotating points round each of the axis
    xaxisrot = np.array([[np.cos(alpha),-np.sin(alpha)],[np.sin(alpha),np.cos(alpha)]])
    yaxisrot = np.array([[np.cos(beta),np.sin(beta)],[-np.sin(beta),np.cos(beta)]])
    #rotate points
    xalpha = np.dot(x,xaxisrot)
    xbeta = np.dot(xalpha,yaxisrot)
    return xbeta 
    
def mate(xa,xb):
    xaxbx = []
    xaxby = []
    for i in range(0,n):
        if xa[i,1] >= 0.0:
            xaxbx.append(x1[i,0])
            xaxby.append(x1[i,1])
        if xb[i,1] < 0.0:
            xaxbx.append(x2[i,0])
            xaxby.append(x2[i,1])      
    xaxb = np.array([xaxbx,xaxby]).transpose(1,0) 
    #make sure there are n points after we have mated
    if len(xaxb) > n:
        lendiff = len(xaxb) - n
        xaxb = xaxb[:-lendiff]
    if len(xaxb) < n:
        lendiff = n - len(xaxb)
        newrows = 2.0*np.random.random((lendiff,2))-1.0 #!!!
        xaxb = np.concatenate((xaxb,newrows))   
        xaxb = proj(xaxb)                               
    finalenergy = calcenergy(xaxb)
    return xaxb, finalenergy

def materelax(xa, xb):
    xaxb, eaeb = mate(xa, xb)
    xaxb, eaeb = relax(xaxb,0.2,eaeb)
 
    return xaxb, eaeb
    

def geneticalgorithm(x1,x2,x3,x4):  
    #rotate spheres by a random amount  
    x1rot = randomrotation(x1)
    x2rot = randomrotation(x2)
    x3rot = randomrotation(x3)
    x4rot = randomrotation(x4) 
    
    #Mate and relax
    x1x2, e1e2 = materelax(x1rot,x2rot)
    x1x3, e1e3 = materelax(x1rot,x3rot)
    x1x4, e1e4 = materelax(x1rot,x4rot)
    x2x1, e2e1 = materelax(x2rot,x1rot)
    x2x3, e2e3 = materelax(x2rot,x3rot)
    x2x4, e2e4 = materelax(x2rot,x4rot)
    x3x1, e3e1 = materelax(x3rot,x1rot)
    x3x2, e3e2 = materelax(x3rot,x2rot)
    x3x4, e3e4 = materelax(x3rot,x4rot)
    x4x1, e4e1 = materelax(x4rot,x1rot)
    x4x2, e4e2 = materelax(x4rot,x2rot)
    x4x3, e4e3 = materelax(x4rot,x3rot)
    
    #find lowest 4 energies
    
    #order energies
    elist = [e1e2,e1e3,e1e4,e2e1,e2e3,e2e4,e3e1,e3e2,e3e4,e4e1,e4e2,e4e3,finalenergy1,finalenergy2,finalenergy3,finalenergy4]
    xlist = [x1x2,x1x3,x1x4,x2x1,x2x3,x2x4,x3x1,x3x2,x3x4,x4x1,x4x2,x4x3,x1,x2,x3,x4]
    newlist = [e1e2]
    
    #make list of all different configurations and positions in our list
    for i in range(0,len(elist)):
        count = 0
        for j in range(0,len(newlist)):
            difff = abs(elist[i] - newlist[j])
            if difff>0.000001: #only compares against one before adding
                count = count +1
                if count == len(newlist):
                    newlist.append(elist[i])
    #print("list")
    #print newlist
                    
    #reorder both lists                 
    newlist, xlist = bubblesort(newlist, xlist)
    
    higheste.append(newlist[-1:])
    
    #if less than 4 minima then duplicate some spheres
    pointi = 2
    while len(newlist) < 4:
        pointi = pointi + 1
        newlist.append(elist[pointi])
        xlist.append(xlist[pointi])
    
    #take lowest 4 energies
    opt4 = newlist[:4]
    points4 = xlist[:4]
   
    x11 = points4[0] #points of lowest energy configuration
    x22 = points4[1]
    x33 = points4[2]
    x44 = points4[3]
    
    e11 = opt4[0] #lowest energy
    e22 = opt4[1] #2nd lowest energy
    e33 = opt4[2] #3rd lowest
    e44 = opt4[3] #4th lowest
    
    loweste.append(e11)
    lowestx.append(x11)
    
    return x11,x22,x33,x44,e11,e22,e33,e44
    
points = open('poadpoints17.txt','w+')
energies = open('poadenergy17.txt','w+')
#points = open('poadpoints30.txt','w+')
#energies = open('poadenergy30.txt','w+')
#points = open('poadpoints38.txt','w+')
#energies = open('poadenergy38.txt','w+')
#points = open('poadpoints56.txt','w+')
#energies = open('poadenergy56.txt','w+')
#points = open('poadpoints69.txt','w+')
#energies = open('poadenergy69.txt','w+')
                
# set the number of points: 2,12,17,30,38,56,69
n=17


looplist=[]
energylist=[]
loweste = []
lowestx = []
higheste = [] 

"""Set up 4 discs"""

#assign random start points on four spheres
random.seed()
x1=initialproj((2.0*np.random.random((n,2))-1.0))
x2=initialproj((2.0*np.random.random((n,2))-1.0))
x3=initialproj((2.0*np.random.random((n,2))-1.0))
x4=initialproj((2.0*np.random.random((n,2))-1.0))
    
# calculate the initial energy
energy1=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance1=np.sqrt(sum((x1[i]-x1[j])**2)) 
        if distance1<>0:
            energy1=energy1+1.0/distance1
energy2=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance2=np.sqrt(sum((x2[i]-x2[j])**2))
        if distance2<>0:
            energy1=energy1+1.0/distance1
        energy2=energy2+1.0/distance2
energy3=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance3=np.sqrt(sum((x3[i]-x3[j])**2))
        if distance3<>0:
            energy1=energy1+1.0/distance1
        energy3=energy3+1.0/distance3
energy4=0.0
for i in range(0,n):
    for j in range(i+1,n):
        distance4=np.sqrt(sum((x4[i]-x4[j])**2))
        if distance4<>0:
            energy1=energy1+1.0/distance1
        energy4=energy4+1.0/distance4

x1, finalenergy1 = relax(x1,0.2,energy1)
x2, finalenergy2 = relax(x2,0.2,energy2)
x3, finalenergy3 = relax(x3,0.2,energy3)
x4, finalenergy4 = relax(x4,0.2,energy4)

# output final energy to the screen and points to a file
#print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy1))
#print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy2))
#print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy3))
#print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy4))

#run through algortihm 5 times
for iterations in range(0,5):
    x1,x2,x3,x4,finalenergy1,finalenergy2,finalenergy3,finalenergy4 = geneticalgorithm(x1,x2,x3,x4)   

empty1=[]
for e in range(0,len(higheste)):
    empty1.append(0)

#output the lowest and highest energy for each value of n
highestenergy, empty1 = bubblesort(higheste, empty1)
lowestenergy, lowestx = bubblesort(loweste, lowestx)

hevalue = highestenergy[-1:]
levalue = lowestenergy[:1]
hh = hevalue[0]

finalenergy2 = levalue[0]

x = lowestx[0]

print("Final energy of {0} ".format(n) + "points = {0:.6f} \n".format(finalenergy2))


for i in range(0,n): 
    points.write('%s %s \n' % (x[i,0],x[i,1]))
energies.write('%s %s %s \n' % (n, levalue[0], hh[0]))

#print("Completed {0} ".format(n) + "points = {0:.6f} \n".format(levalue[0]))
    
endtime = time.time()
#print(endtime-starttime)


points.close()
energies.close()  


