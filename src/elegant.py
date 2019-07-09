
'''
Created on Feb 12, 2019

@author: rravell
'''
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt

if __name__ == '__main__':
    
    outputsAlice = [3,3,3,5]
    outputsBob = [3,3,3,3]
    
    
    psi=createMaxEntState(2)
    
    aliceEffects=[[],[],[],[]]
    aliceEffects[0]=[1/2*(qt.qeye(2)+qt.sigmaz()),1/2*(qt.qeye(2)-qt.sigmaz()),0*qt.qeye(2)]
    aliceEffects[1]=[1/2*(qt.qeye(2)+qt.sigmax()),1/2*(qt.qeye(2)-qt.sigmax()),0*qt.qeye(2)]
    aliceEffects[2]=[1/2*(qt.qeye(2)+qt.sigmay()),1/2*(qt.qeye(2)-qt.sigmay()),0*qt.qeye(2)]
    aliceEffects[3]=[1/4*(qt.qeye(2)-1/np.sqrt(3)*(qt.sigmaz()+qt.sigmax()+qt.sigmay())),1/4*(qt.qeye(2)-1/np.sqrt(3)*(qt.sigmaz()-qt.sigmax()-qt.sigmay())),1/4*(qt.qeye(2)+1/np.sqrt(3)*(qt.sigmaz()-qt.sigmax()+qt.sigmay())),1/4*(qt.qeye(2)+1/np.sqrt(3)*(qt.sigmaz()+qt.sigmax()-qt.sigmay())),0*qt.qeye(2)]
    
    bobEffects=[[],[],[],[]]
    bobEffects[0]=[1/2*(qt.qeye(2)+1/np.sqrt(3)*(qt.sigmaz()+qt.sigmax()-qt.sigmay())),1/2*(qt.qeye(2)-1/np.sqrt(3)*(qt.sigmaz()+qt.sigmax()-qt.sigmay())),0*qt.qeye(2)]
    bobEffects[1]=[1/2*(qt.qeye(2)+1/np.sqrt(3)*(qt.sigmaz()-qt.sigmax()+qt.sigmay())),1/2*(qt.qeye(2)-1/np.sqrt(3)*(qt.sigmaz()-qt.sigmax()+qt.sigmay())),0*qt.qeye(2)]
    bobEffects[2]=[1/2*(qt.qeye(2)+1/np.sqrt(3)*(-qt.sigmaz()+qt.sigmax()+qt.sigmay())),1/2*(qt.qeye(2)-1/np.sqrt(3)*(-qt.sigmaz()+qt.sigmax()+qt.sigmay())),0*qt.qeye(2)]
    bobEffects[3]=[1/2*(qt.qeye(2)+1/np.sqrt(3)*(-qt.sigmaz()-qt.sigmax()-qt.sigmay())),1/2*(qt.qeye(2)-1/np.sqrt(3)*(-qt.sigmaz()-qt.sigmax()-qt.sigmay())),0*qt.qeye(2)]
    #bobEffects[4]=[1/2*(qt.qeye(2)+1/np.sqrt(2)*(qt.sigmay()+qt.sigmaz())),1/2*(qt.qeye(2)-1/np.sqrt(2)*(qt.sigmay()+qt.sigmaz())),0*qt.qeye(2)]
    #bobEffects[5]=[1/2*(qt.qeye(2)+1/np.sqrt(2)*(qt.sigmay()-qt.sigmaz())),1/2*(qt.qeye(2)-1/np.sqrt(2)*(qt.sigmay()-qt.sigmaz())),0*qt.qeye(2)]
    
    dist=[]
    for x in range (0,len(outputsAlice)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAlice[x]):
                for b in range (0,outputsBob[y]):
                    dist.append((qt.tensor(aliceEffects[x][a],bobEffects[y][b])*psi*psi.dag()).tr())
     
     
#     vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #NO DISYUNTOS
#     vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #SI DISYUNTOS
#     Disjoint=np.concatenate((vector2,vector2,vector2,vector2,vector2,vector1,vector2,vector1,vector2,vector2,vector1,vector1,vector2,vector1,vector1,vector1))
#     
#     dist=Disjoint
#     
    matrixHRepresentation = createConstraintsMatrixForEfficiencyDual(outputsAlice, outputsBob) 
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = tuple(np.concatenate(([0],dist)))
    
    lp = cdd.LinProg(mat)
    lp.solve()
    print(lp.status)
    print(lp.obj_value)
    print(" ".join("{0}".format(val) for val in lp.primal_solution))
    
        
    
