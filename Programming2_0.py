# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 12:16:10 2020

@author: Owner
"""
import numpy as np
import matplotlib.pyplot as plt
import re, os 
import math
from scipy.optimize import curve_fit

#creating intial lists to append
r_=[]
theta_=[]
e_=[]
r_1=[]
theta_1=[]
E_1=[]

#User input to direct to the correct folder
nmdir=str(input("Directory: "))

#loop through every file in folder
#findall recognises certain pattern of numbers identifying r and theta from the filename
for filename in os.listdir(nmdir):
    r_t_e=re.findall(r"\d*\.\d+|\d+",str(filename))
    r_.append(float(r_t_e[1]))
    theta_.append(float(r_t_e[2]))
              
    #open each file then search for the line containing energy to then extract the energy
    f=open(nmdir + '/'+str(filename),"r")
    for line in f:
        if "SCF Done:" in line:
            l=line.split()
            e_.append(float(l[4]))
    f.close()

        
#converting lists into numpy arrays     
r=np.array(r_)
theta=np.array(theta_)
E=np.array(e_)  

#Plotting
fig=plt.figure()
ax=plt.axes(projection='3d')

#plotting the energy surface
ax.plot_trisurf(r,theta,E,cmap='jet')

#customisation
ax.set_title('Energy surface for geometry (r,θ)')
ax.set_xlabel('Bond length/ Å')
ax.set_ylabel('θ/ °')
ax.set_zlabel('Energy/ Hartrees')


#setting X as a 2D array to be used in optimisation and curve fitting
X=r,theta

#Run curve fitting once
#approximating this curve with equation given in notes 
def objective(X,E0,k_r,r_0,k_theta,theta_0):
    return E0+0.5*k_r*(r-r_0)**2+0.5*k_theta*(theta-theta_0)**2
fit_1,_=curve_fit(objective,X,E) 

#fit_1 contains all the constants from the approximation

#2nd run through of fit
#Expand approximation around minima determined to get a better approximation of the turning point

for i in range(len(r)):
    if r[i] > fit_1[2]-0.15 and r[i] < fit_1[2]+0.15:           #restricting the values of r to between 0.15 of the minima
        if theta[i] > fit_1[4]-15 and theta[i] < fit_1[4]+15:       #restrict values of theta to between 15 of minima
            r_1.append(r[i])        #create new lists to run second approximation
            theta_1.append(theta[i])
            E_1.append(E[i])
#define new 2D array to be used in optimisation and curve fitting
X_1=(r_1,theta_1)
#using equation found in notes for curve fitting
def objective(X_1,E0,k_r,r_0,k_theta,theta_0):
    return E0+0.5*k_r*(r_1-r_0)**2+0.5*k_theta*(theta_1-theta_0)**2
fit_2,_=curve_fit(objective,X_1,E_1) 
#fit_2 returns all values of the equation 
print(str(fit_2))
K_1=fit_2[0]+0.5*fit_2[1]*(r_1-fit_2[2])**2+0.5*fit_2[3]*(theta_1-fit_2[4])**2

#calculation of frequencies
#constants
mu=1.66*(10**-27)
c=2.998*(10**10)
a=4.35974*10**(-18)     #conversion of hartree to joules
b=1/10**(-20)       #conversion from (Angstrom)^-2 to m^-2

#equation for frequencies from the notes
v1=(1/(2*math.pi))*((fit_2[1]*a*b/(2*mu))**(1/2))/c
v2=(1/(2*math.pi))*(180/math.pi)*((fit_2[3]*a/(0.5*mu*((fit_2[2]*10**(-10))**2)))**(1/2))/c
print('Vibrational frequencies')
print('Symmetric stretch: '+str(round(v1,2)) +' cm-1')
print('Bend: '+str(round(v2,3))+' cm-1')


