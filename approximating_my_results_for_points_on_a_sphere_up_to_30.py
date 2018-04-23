# graph of my found energy values against N, and an approximation of the curve
# graph showing the differeces between my results and the approximation
# these are coloured accoring to symmetries

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

data = np.genfromtxt('energies30.txt')

nopf, energyvalf = np.hsplit(data,2)
nop = nopf.flatten()
energyval = energyvalf.flatten()

def f(n,alpha,beta):
    return 0.5*n**2.0 - alpha*n**1.5 + beta*n**0.5

guess = np.array([0.0,0.0])
optimise = opt.curve_fit(f,nop,energyval,p0=guess)[0]
a = optimise[0]
b = optimise[1]

fit = f(nop,a,b)

residual = energyval - fit

#plot y=0 line
zerolinex = np.linspace(0,35,num = 2)
zeroliney = [0,0]



#colouring odd and even
resodd = []
reseven = []
nopodd = []
nopeven = []
for i in range(0,1+len(residual)/2):
    reseven.append(residual[2*i])
    nopeven.append(nop[2*i])
    resodd.append(residual[(2*i)-1])
    nopodd.append(nop[(2*i)-1])

#colouring symmetries
nopp = [4,6,12,16,22,24,28] 
nopc = [11,13,19,21,25,26]
nopd = [2,3,5,7,8,9,10,14,15,17,18,20,23,27,29,30]
resp=[]
resd=[]
resc=[]

for i in range(0,len(nopp)):
    for j in range(0,len(nop)):
        if nopp[i]==nop[j]:
            resp.append(residual[j])
for i in range(0,len(nopd)):
    for j in range(0,len(nop)):
        if nopd[i]==nop[j]:
            resd.append(residual[j])
for i in range(0,len(nopc)):
    for j in range(0,len(nop)):
        if nopc[i]==nop[j]:
            resc.append(residual[j])
            
    





fig1 = plt.figure(1)
fig2 = plt.figure(2)

ax5=fig1.add_subplot(111)
ax5.plot(nop, fit, linestyle='-',color = 'magenta', linewidth=1.0, label='Approximation')
ax5.plot(nop, energyval, linestyle='', marker='x', c='teal', mew=1.5, label='My results')
ax5.legend(loc='upper left', fontsize='12')
ax5.set_xlabel("Number of points")
ax5.set_ylabel("Energy")


ax6=fig2.add_subplot(212)
plt.plot(zerolinex, zeroliney, linestyle='-', c = 'b')
ax7=fig2.add_subplot(211)
plt.plot(zerolinex, zeroliney, linestyle='-', c = 'b')

ax7.plot(nopodd, resodd, linestyle='', marker='o', c='lightblue', label='odd N')
ax7.plot(nopeven, reseven, linestyle='', marker='o', c='blue', label='even N')
ax6.plot(nopp, resp, linestyle='', marker='o', c='darkgreen', label='Platonic')
ax6.plot(nopd, resd, linestyle='', marker='o', c='mediumspringgreen', label='Dihedral')
ax6.plot(nopc, resc, linestyle='', marker='o', c='azure', label='Cylic')

ax6.set_xlim([0,32])
ax6.set_ylim([-0.11,0.09])
ax6.legend(loc='upper left', fontsize='8')
ax6.set_xlabel("Number of points")
ax6.set_ylabel("Actual energy - approximated energy")


ax7.set_xlim([0,32])
ax7.set_ylim([-0.11,0.09])
ax7.legend(loc='upper left', fontsize='8')
ax7.set_xlabel("Number of points")
ax7.set_ylabel("Actual energy - approximated energy")



plt.show() 

