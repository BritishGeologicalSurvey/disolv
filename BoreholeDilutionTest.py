# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:22:42 2018

@author: sarcol
"""

import os
import numpy as np
import SolveEquation
import matplotlib.pyplot as plt

def ConvertFEC(InData,Temp):
    
    FEC20 = (1 + 0.024 * (Temp - 20)) * InData
    Concentrations = (1870 - np.sqrt(1870**2 - 80 * FEC20))/80

    return Concentrations    


#---------------------Get input parameters-----------------------

InDat = os.path.join("Input","in.dat")

Parameters = np.genfromtxt(InDat,skip_footer=1)

length = Parameters[0]
z = Parameters[1]
A = Parameters[2]
alpha = Parameters[3]
Cc = Parameters[4]
FECconvert = Parameters[5]
Temp = Parameters[6]
Calibrate = Parameters[7]

Bounds = np.genfromtxt(InDat,skip_header=16,skip_footer=1)
t = np.genfromtxt(InDat,skip_header=18)
t = np.concatenate((np.array([0]),t))

#----------------------------Grid---------------------------------

N_nodes = int(round((length+z)/z))
x = np.linspace(0,length,N_nodes)

#---------Initial condition and Observations-----------------------

in_con_raw = np.genfromtxt(os.path.join("Input","initial_condition.dat"),skip_header=1)
ObsProfilesRaw = np.genfromtxt(os.path.join("Input","measuredprofiles.dat"),skip_header=1)

x = np.linspace(0,length,N_nodes)

NObs = int(np.shape(ObsProfilesRaw)[1]/2)
ObservedProfiles = np.zeros([len(x),NObs])

if Parameters[5] == 2:
    Concentrations = ConvertFEC(in_con_raw[:,1],Temp)
    for i in range(NObs):
        ObsProfilesRaw[:,2*i + 1] = ConvertFEC(ObsProfilesRaw[:,2*i + 1],Temp)
else:
    Concentrations = in_con_raw[:,1]

in_con = np.interp(x,in_con_raw[:,0],Concentrations)

#Oberservation x points for automatic calibration

for i in range(NObs):
    ObservedProfiles[:,i] = np.interp(x,ObsProfilesRaw[:,2*i],ObsProfilesRaw[:,2*i + 1])



#------------------------Make sure Qin = Qout---------------------

indata = np.genfromtxt(os.path.join("Input","flows.dat"),skip_header=1)[:,:2]

Nflows = np.shape(indata)[0]

for i in range(Nflows):
    
    indata[i,0] = round(indata[i,0],1)

inflows = np.zeros([Nflows,2])
outflows = np.zeros([Nflows,2])

for i in range(Nflows):
    
    if indata[i,1] < 0:
        
        outflows[i,0] = indata[i,0]
        outflows[i,1] = -indata[i,1]
    
    if indata[i,1] > 0:
        
        inflows[i,0] = indata[i,0]
        inflows[i,1] = indata[i,1]
        
 
#-------------------Call function that solves equation----------------
if Calibrate == 0:
    sim_profiles = SolveEquation.forward(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A)


#-------------------------Automatic calibration------------------------
 
if Calibrate == 1:
    
    output = SolveEquation.inverse(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A,ObservedProfiles,Bounds)
    
    alpha = output.x[0]
    inflows[:,1] = output.x[1:(Nflows+1)]/np.sum(output.x[1:(Nflows+1)]) * output.x[-1]
    outflows[:,1] = output.x[(Nflows+1):-1]/np.sum(output.x[(Nflows+1):-1])* output.x[-1]
    
    sim_profiles = SolveEquation.forward(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A)
    
             
    
out_values = sim_profiles[1:,:]
OutStr = []
for i in range(len(t)-1):
    OutStr.append(str(t[i+1])+'\t')
OutStr.append('\n')
for j in range(N_nodes):
    for i in range(len(t)-1):
        OutStr.append(str(out_values[i,j]) + '\t')
    OutStr.append('\n')

f = open(os.path.join("Output","profiles.out"), 'w')
for item in OutStr:
	f.write("%s" % item)
f.close()

#-----------------Plot results------------------

plt.figure()
for i in range(NObs):
    plt.plot(sim_profiles[i+1,:],x,label=str(t[i+1]))
    plt.scatter(ObsProfilesRaw[:,2*i + 1],ObsProfilesRaw[:,2*i])

plt.legend()
plt.xlabel('Salinity (kg/m$^3$)')
plt.ylabel('Depth below ground (m)')
plt.gca().invert_yaxis()
plt.savefig(os.path.join('Output', 'result.png'))