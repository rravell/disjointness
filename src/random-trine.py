
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

def createRandomTrine():
    phase=np.random.uniform(-np.pi/2,np.pi/2)
    vectors=[[np.sin(phase+k*2*np.pi/3),0,np.cos(phase+k*2*np.pi/3)] for k in range(0,3)]
    return list(map(lambda vector: effectForQubitPovm(1/3, createQubitObservable(vector)),vectors))

if __name__ == '__main__':
    
    outputsAlice = [4,4,4]
    outputsBob = [4,4,4,4]
    
    aliceTrines = [createRandomTrine() for _ in range(0,len(outputsAlice))]
    bobTrines = [createRandomTrine() for _ in range(0,len(outputsBob))]
        
    aliceEffects = addEffectsForAbortOutcomes(aliceTrines,2)
    bobEffects = addEffectsForAbortOutcomes(bobTrines,2)
    
    psi=createMaxEntState(2)
    
    dist=computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)

    matrixHRepresentation = createConstraintsMatrixForEfficiencyDual(outputsAlice, outputsBob) 
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = tuple(np.concatenate(([0],dist)))
    
    lp = cdd.LinProg(mat)
    lp.solve()
    print(lp.status)
    print(lp.obj_value)
    print(" ".join("{0}".format(val) for val in lp.primal_solution))
    
        
    
