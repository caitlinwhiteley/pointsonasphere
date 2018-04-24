# python code to compute points in a disk
# uses method of steepest decent
# draws the points on the disk 
# imports my results for N=17,30,38,56,69, which I found using the genetic algorithm


import numpy as np
import matplotlib.pyplot as plt
import time
from sys import exit
 

starttime = time.time()


"""FUNCTIONS"""
def proj(x,n) :
    for i in range(0,n):
        norm=np.sqrt(sum(x[i]**2))
        if norm>1:
            x[i]=x[i]/norm
    return x
    
def calcenergy(x,n):
    energyval=0.0
    for i in range(0,n):
        for j in range(i+1,n):
            distance=np.sqrt(sum((x[i]-x[j])**2))
            if distance > 0:
                energyval=energyval+1.0/distance
    return energyval
    
def steepestdecent(x, amplitude,n): #called within relax function
    forceslist = [] #empty list we will add the resultant components of each xi to
    forces_pointsjoni = np.zeros([n,2]) #empty array which we will add the resultant force of each point on a single other point
    for i in range(0,n):
        for j in range(0,n):
            if x[i,0]-x[j,0]==0 and x[i,1]-x[j,1]==0 and x[i,1]-x[j,1]==0:
                forces_pointsjoni[j] = 0
            else:
                forces_pointsjoni[j] = ((x[i]-x[j])/(np.sqrt(sum((x[i]-x[j])**2))**3)) #resultant forces of point j on point i
        resultanttotals = np.nansum(forces_pointsjoni, axis=0) #sum up the columns of get total x,y,z compnents of forces acting on i
        forcetilda = resultanttotals - np.dot(resultanttotals, x[i]) * x[i] #taking dot product to only take forces in relevant direction
        forceslist.append(forcetilda) #adding our forcetilda to our empty list
    for i in range(0,n):
        x[i] = x[i] + amplitude*forceslist[i]
    return x
    
def relax(x,amplitude,energy,n):
    loop=0
    while amplitude > 0.00005: 
        loop=loop+1 
        #store the value of the energy to compare to when we calucalte a new energy
        oldenergy = energy.copy()
        # store the old coordinates of this point
        oldx = x.copy()
        #calculate resultant forces of each point
        x = steepestdecent(x, amplitude,n)
        x = proj(x,n)
        # calculate the energy
        energy = calcenergy(x,n)        
        #reducing step size if energy increases 
        difference = energy-oldenergy
        if(difference>0.0): #if the old energy is lower than the new energy reduce step size
            amplitude = amplitude/2.0
            x = oldx.copy()
            energy = oldenergy.copy()      
    return x, energy

def bubblesort1(newlist):
    for ii in range(len(newlist)-1,0,-1):
        for i in range(ii):
            if newlist[i]>newlist[i+1]:
                temp = newlist[i]
                newlist[i] = newlist[i+1]
                newlist[i+1] = temp
    return newlist

def bubblesort2(newlist,xlist):
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
    

    


"""INITIAL SETUP"""                                                
                                                                                
# set the number of points
n = 70



#Find an initial configuration, likely to be incorrect but we'll use it later    
looplist=[]
energylist=[]


"""N<12: all points on boundary, 11<n<17 all points on boundary except one
in the middle"""

#If 1<n<12, all points on boundary
if 2 <= n <= 16:
    if n<=11:        
        x = np.zeros((n,2))
        for b in range(0,n):
            phi0 = 2.0*b*np.pi/(n)
            x[b] = (np.cos(phi0),np.sin(phi0))   
    else:
        x = np.zeros((n,2))
        for b in range(0,n-1):
            phi0 = 2.0*b*np.pi/(n-1)
            x[b] = (np.cos(phi0),np.sin(phi0)) 
    
    energy = calcenergy(x,n)  

    # output final energy to the screen and points to a file
    print("Energy of {0} ".format(n) + "points = {0:.6f} \n".format(energy))
    
    #PLOT
    #Create a disk
    theta = np.linspace(0,2.0*np.pi,400)
    xcoord = np.cos(theta)
    ycoord = np.sin(theta)
    #convert data
    x1=[]
    x2=[]
    for i in range(0,n):
        x1.append(x[i,0])
        x2.append(x[i,1])    
    #Render
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.scatter(xcoord,ycoord, color='blue', s=1)
    ax.scatter(x1,x2,s=80, color='black')
    plt.axis('off')
    ax.set_aspect("equal")
    ax.set_title("{0} ".format(n)+"points on a disc")
    plt.show()
    
    exit()


"""Terminate if x is out of bounds"""

#if n=1 or n>89, terminate program
if n>89:
    print 'Cannot compute, choose 1<N<90'
    exit()
if n==1:
    print 'Cannot compute, choose 1<N<90'
    exit()
    
    
"""For n within the range of my program, import points of the smallest n 
within the group of concentric rings. Previously found using the genetic
algorithm and saved to a text file"""   
#for 16<n<90, compute
    
if 17 <= n <= 29:
    nmin = 17
    data = np.genfromtxt('poadpoints17.txt')
if 30 <= n <= 37:
    nmin = 30
    data = np.genfromtxt('poadpoints30.txt')   
if 38 <= n <= 55:
    nmin = 38
    data = np.genfromtxt('poadpoints38.txt')    
if 56 <= n <= 68:
    nmin = 56
    data = np.genfromtxt('poadpoints56.txt')    
if 69 <= n <= 89:
    nmin = 69
    data = np.genfromtxt('poadpoints69.txt')      

x0, x1 = np.hsplit(data,2)
x0 = x0.flatten()
x1 = x1.flatten()

#create an array of the points of the minimum n in the group
xmin = np.zeros((len(x0),2))
for i in range(0,len(x0)):
    xmin[i,0] = x0[i]
    xmin[i,1] = x1[i]

xk = xmin.copy()



for k in range(0,n-nmin+1):

    #calculate radius of each ring of previous configurations, add to list
    ringslist = [1]
    xkr = []
    
    ringslist = [1]
    xkr = []
    for i in range(0,nmin+k-1):
        ringr = np.sqrt(sum((xk[i])**2))
        xkr.append(ringr)
        if ringr > 0.005:
            count = 0
            for j in range(0,len(ringslist)):
                if abs(ringslist[j] - ringr) > 0.2:
                    count = count + 1
                if count == len(ringslist):
                    ringslist.append(ringr)
            
  
    z = np.zeros((1,len(ringslist)))
    for i in range(0,len(ringslist)):
        for j in range(0,len(xkr)):
           if abs(xkr[j] - ringslist[i]) < 0.1:
                z[0,i] = z[0,i] + 1
    listt=[]       
    for i in range(0,len(ringslist)):
        listt.append(int(z[0,i]))

    

        
    listt = bubblesort1(listt)
    ringslist = bubblesort1(ringslist)
    
    
    
    #add a point to each radius in turn
    if len(listt)==2:
        list1 = [listt[0]+1,listt[1]]
        list2 = [listt[0],listt[1]+1]
    if len(listt)==3:
        list1 = [listt[0]+1,listt[1],listt[2]]
        list2 = [listt[0],listt[1]+1,listt[2]]
        list3 = [listt[0],listt[1],listt[2]+1]
    if len(listt)==4:
        list1 = [listt[0]+1,listt[1],listt[2],listt[3]]
        list2 = [listt[0],listt[1]+1,listt[2],listt[3]]
        list3 = [listt[0],listt[1],listt[2]+1,listt[3]]
        list4 = [listt[0],listt[1],listt[2],listt[3]+1]
    
    xlist = []
    elist = []

    
    x = np.zeros((nmin+k,2))
    for b in range(0,list1[0]):
        r = ringslist[0]
        phi0 = 2.0*b*np.pi/(list1[0])
        x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    for b in range(list1[0],list1[1]+list1[0]):
        r = ringslist[1]
        phi0 = 2.0*b*np.pi/(list1[1])
        x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    if len(listt) > 2.5:
        for b in range(list1[0]+list1[1],list1[0]+list1[1]+list1[2]):    
            r = ringslist[2]
            phi0 = 2.0*b*np.pi/(list1[2])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    if len(listt) == 4:
        for b in range(list1[0]+list1[1]+list1[2],list1[0]+list1[1]+list1[2]+list1[3]):    
            r = ringslist[3]
            phi0 = 2.0*b*np.pi/(list1[3])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    energy1 = calcenergy(x,nmin+k)     
    x,energy = relax(x,0.05,energy1,nmin+k)
    
    
    xlist.append(x)
    elist.append(energy)

    x = np.zeros((nmin+k,2))
    for b in range(0,list2[0]):
        r = ringslist[0]
        phi0 = 2.0*b*np.pi/(list2[0])
        x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    for b in range(list2[0],list2[1]+list2[0]):
        r = ringslist[1]
        phi0 = 2.0*b*np.pi/(list2[1])
        x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    if len(listt) > 2.5:
        for b in range(list2[0]+list2[1],list2[0]+list2[1]+list2[2]):    
            r = ringslist[2]
            phi0 = 2.0*b*np.pi/(list2[2])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    if len(listt) == 4:
        for b in range(list2[0]+list2[1]+list2[2],list2[0]+list2[1]+list2[2]+list2[3]):    
            r = ringslist[3]
            phi0 = 2.0*b*np.pi/(list2[3])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
    energy1 = calcenergy(x,nmin+k)
    x,energy = relax(x,0.05,energy1,nmin+k)
    
    xlist.append(x)
    elist.append(energy)
    
    if len(listt) > 2.5:
        x = np.zeros((nmin+k,2))
        for b in range(0,list3[0]):
            r = ringslist[0]
            phi0 = 2.0*b*np.pi/(list3[0])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        for b in range(list3[0],list3[1]+list3[0]):
            r = ringslist[1]
            phi0 = 2.0*b*np.pi/(list3[1])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        for b in range(list3[0]+list3[1],list3[0]+list3[1]+list3[2]):
            r = ringslist[2]
            phi0 = 2.0*b*np.pi/(list3[2])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        if len(listt) == 4:
            for b in range(list3[0]+list3[1]+list3[2],list3[0]+list3[1]+list3[2]+list3[3]):
                r = ringslist[3]
                phi0 = 2.0*b*np.pi/(list3[3])
                x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        energy1 = calcenergy(x,nmin+k)
        x,energy = relax(x,0.05,energy1,nmin+k)
        
        xlist.append(x)
        elist.append(energy)

    if len(listt) == 4:
        x = np.zeros((nmin+k,2))
        for b in range(0,list4[0]):
            r = ringslist[0]
            phi0 = 2.0*b*np.pi/(list4[0])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        for b in range(list4[0],list4[1]+list4[0]):
            r = ringslist[1]
            phi0 = 2.0*b*np.pi/(list4[1])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        for b in range(list4[0]+list4[1],list4[0]+list4[1]+list4[2]):
            r = ringslist[2]
            phi0 = 2.0*b*np.pi/(list4[2])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        for b in range(list4[0]+list4[1]+list4[2],list4[0]+list4[1]+list4[2]+list4[3]):
            r = ringslist[3]
            phi0 = 2.0*b*np.pi/(list4[3])
            x[b] = (r*np.cos(phi0),r*np.sin(phi0))
        energy1 = calcenergy(x,nmin+k)
        x,energy = relax(x,0.05,energy1,nmin+k)
        
        xlist.append(x)
        elist.append(energy)
    
    

    mine, minx = bubblesort2(elist,xlist)
    minx=minx[0]

    
    xk = minx.copy()
    
    nk = nmin+k
    print("Energy of {0} ".format(nk) + "points = {0:.6f} \n".format(mine[0]))



#PLOT
#Create a disk
theta = np.linspace(0,2.0*np.pi,400)
xcoord = np.cos(theta)
ycoord = np.sin(theta)
#convert data
x1=[]
x2=[]
for i in range(0,nmin+k):
    x1.append(xk[i,0])
    x2.append(xk[i,1])    
#Render
fig = plt.figure()
ax = fig.add_subplot(211)
ax.scatter(xcoord,ycoord, color='blue', s=1)
ax.scatter(x1,x2,s=80, color='black')
plt.axis('off')
ax.set_aspect("equal")
ax.set_title("{0} ".format(nmin+k)+"points on a disc")
plt.show()

endtime = time.time()
#print(endtime-starttime)


