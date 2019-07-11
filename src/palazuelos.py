
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

def createEffects(observables):
    return list(map(lambda observable :  
                       list(map(lambda eigvector : eigvector*eigvector.dag(),observable)),observables))

def addEffectsForAbortOutcomes(effects):
    return list(map(lambda povm : povm+[0*qt.qeye(dim)],effects))
if __name__ == '__main__':
    
    N=3
    dim=N+1
    outputsAlice = [dim+1 for i in range(0,N)]
    outputsBob = [dim+1 for i in range(0,N)]
    
    
    eigVectAlice= [[] for x in range(1,N+1)]
    
    eigVectBob= [[1/np.sqrt(dim)*np.sum([np.exp(-2*np.pi*1j/dim*q*(r-l/N))*qt.basis(dim, q) for q in range(0,dim)])
                    for r in range(0,dim)] 
                    for l in range(1,N+1)]
    
    aliceEffects = addEffectsForAbortOutcomes(createEffects(eigVectAlice))
    bobEffects = addEffectsForAbortOutcomes(createEffects(eigVectBob))
    psi=createMaxEntState(dim)
     
    distribution = computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
      
    matrixHRepresentation = createConstraintsMatrixForEfficiencyDual(outputsAlice, outputsBob) 
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = tuple(np.concatenate(([0],distribution)))
    
    lp = cdd.LinProg(mat)
    lp.solve()
    print(lp.status)
    print(lp.obj_value)
    print(" ".join("{0}".format(val) for val in lp.primal_solution))
    
        
    
