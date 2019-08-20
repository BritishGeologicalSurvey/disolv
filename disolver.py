# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:22:42 2018

@author: sarcol
"""

import os
import numpy as np
import SolveEquation
import matplotlib.pyplot as plt
import pandas as pd
import argparse



def run(InDir,OutDir, calibrate='False', convertFEC = 'False', method='SLSQP'):
    
    def ConvertFEC(InData,Temp):
    
        FEC20 = (1 + 0.024 * (Temp - 20)) * InData
        Concentrations = (1870 - np.sqrt(1870**2 - 80 * FEC20))/80
    
        return Concentrations    



    #---------------------Get input parameters-----------------------
    
    InDat = os.path.join(InDir,"in.csv")
    
    Parameters = np.genfromtxt(InDat,delimiter=',',skip_footer=1)[:-1,0]
    
    GWlevel = Parameters[0]
    BHdepth = Parameters[1]
    z = Parameters[2]
    A = Parameters[3]
    alpha = Parameters[4]
    Cc = Parameters[5]
    Temp = Parameters[6]
    
    SatColumn = z * round((BHdepth - GWlevel)/z)
    
    
    Bounds = np.genfromtxt(InDat,delimiter=',',skip_header=14,skip_footer=1)[:4]
    t = np.genfromtxt(InDat,delimiter=',',skip_header=16)
    t = t[np.isfinite(t)]
    t = np.concatenate((np.array([0]),t))
    
    #----------------------------Check inputs------------------------
    
    in_con_raw = np.genfromtxt(os.path.join(InDir,"initialcondition.csv"),delimiter=',',skip_header=1)
    
    #---Check for observation data---
    
    ObsFile = os.path.join(InDir,"measuredprofiles.csv")
    ObsExist = os.path.isfile(ObsFile)
    if ObsExist:
        df = pd.read_csv(os.path.join(ObsFile))
        ObsProfilesRaw = df.values
        NObs = int(np.shape(ObsProfilesRaw)[1]/2)
    
    else:
        ObsProfilesRaw = np.copy(in_con_raw)
        ObsProfilesRaw[:,1] = np.NaN
        calibrate='False'
        NObs = 0
    
    
    if len(t)-1 != NObs:
        raise Exception('{} output times given, but '.format(len(t)-1) + str(NObs) + ' observation profiles found. These should be equal.')
    if len(np.shape(in_con_raw)) == 1:
        raise Exception('initialcondition.csv must have two columns: depth and concentration')
    
    
    #---Print information
    if calibrate == 'False':
        print("No automatic calibration")
    elif calibrate == 'True':
        print("Automatic calibration")
    
    print("Equation will be solved at times " + str(t[1:].tolist()).strip('[]'))
    print(str(NObs) + " measured profiles have been found")
    #----------------------------Grid---------------------------------
    
    N_nodes = int((SatColumn+z)/z)
    x = np.linspace(0,SatColumn,N_nodes)
    
    #---------Initial condition and Observations-----------------------
    
    
    ObservedProfiles = np.zeros([len(x),NObs])
    
    if convertFEC == 'True':
        Concentrations = ConvertFEC(in_con_raw[:,1],Temp)
        for i in range(NObs):
            ObsProfilesRaw[:,2*i + 1] = ConvertFEC(ObsProfilesRaw[:,2*i + 1],Temp)
    else:
        Concentrations = in_con_raw[:,1]
    
    in_con = np.interp(x,in_con_raw[:,0]-GWlevel,Concentrations)
    
    #Oberservation x points for automatic calibration
    
    for i in range(NObs):
        ObservedProfiles[:,i] = np.interp(x,ObsProfilesRaw[:,2*i]-GWlevel,ObsProfilesRaw[:,2*i + 1])
    
    
    
    #------------------------Make sure Qin = Qout---------------------
    
    indata = np.genfromtxt(os.path.join(InDir,"flows.csv"),delimiter=',',skip_header=1)
    Nflows = np.shape(indata)[0]
    
    indata[:,0] = indata[:,0] - GWlevel
    indata[:,2:] = indata[:,2:] - GWlevel
    
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
    
    if calibrate == 'False':
        sim_profiles = SolveEquation.forward(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A)
        
    
    #-------------------------Automatic calibration------------------------
     
    elif calibrate == 'True':
        
        FracBounds = ()
        
        if np.shape(indata)[1] == 4:
            for i in range(Nflows):
                FracBounds = FracBounds + ((indata[i,2],indata[i,3]),)
        else:
            for i in range(Nflows):
                FracBounds = FracBounds + ((indata[i,0],indata[i,0]),)
                
        output = SolveEquation.inverse(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A,ObservedProfiles,Bounds,FracBounds,method)
        
        alpha = output.x[0]
        inflows[:,1] = output.x[1:(Nflows+1)]/np.sum(output.x[1:(Nflows+1)]) * output.x[-1]
        outflows[:,1] = output.x[(Nflows+1):-(Nflows+1)]/np.sum(output.x[(Nflows+1):-(Nflows+1)])* output.x[-1]
        inflows[:,0] = output.x[-(Nflows+1):-1]
        outflows[:,0] = output.x[-(Nflows+1):-1]    
        
        sim_profiles = SolveEquation.forward(N_nodes,inflows,outflows,z,alpha,Cc,Nflows,in_con,t,A)
        
        
                 
        
    out_values = np.zeros([N_nodes,NObs])
    
    for j in range(N_nodes):
        for i in range(len(t)-1):
            
            out_values[j,i] = sim_profiles[i+1,j]
    colheads = []
    for i in range(len(t)-1):
        colheads.append("t = " + str(t[i+1]))
    
    dfout = pd.DataFrame(data=out_values,index=x+GWlevel, columns=colheads)
    dfout.index.name = "Depth [L]"
    dfout.to_csv(os.path.join(OutDir,"profiles.csv"))
    
    if calibrate == 'True':
        
        out = "Dispersivity, " + str(alpha) + "\n" + "Flow rates" + "\n" + "Depth [L]" + "," + "Flow [L^2T^-1]" + "\n"
        for i in range(Nflows):
            if inflows[i,1] != 0:  
                out = out + str(inflows[i,0] + GWlevel) + ',' + str(inflows[i,1]) + '\n' 
            else:
                out = out + str(outflows[i,0] + GWlevel) + ',' + str(outflows[i,1]) + '\n'
        f = open(os.path.join(OutDir,'Output.csv'), 'w')
        for item in out:
        	f.write("%s" % item)
        f.close()    
                
    
    #-----------------Plot results------------------
    
    plt.figure(figsize=(5,8))
    for i in range(len(t)-1):
        plt.plot(sim_profiles[i+1,:],x+GWlevel,label=str(t[i+1]))
        if ObsExist:
            plt.scatter(ObsProfilesRaw[:,2*i + 1],ObsProfilesRaw[:,2*i])
    
    plt.legend()
    plt.xlabel('Salinity (kg/m$^3$)')
    plt.ylabel('Depth below ground (m)')
    plt.gca().invert_yaxis()
    plt.savefig(os.path.join(OutDir, 'profiles.png'))


def main():
    """
    Run the solver with default arguments
    """
    in_dir = os.path.join(os.getcwd(), 'Input')
    out_dir = os.path.join(os.getcwd(), 'Output')
    run(in_dir, out_dir, calibrate='False', convertFEC = 'False', method='SLSQP')
    
    
    
if __name__ == "__main__":
    
    main()
    
