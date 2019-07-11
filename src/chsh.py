
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
    
    outputsAlice = [3,3,5]
    outputsBob = [3,3,4,5]
    
    aliceUnBlochVectors = [[1,0,0],[0,1,0]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceUnBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
    
    delta=0.01
    tetrahedronAlice = [[delta,0,1],[-delta,0,1],[0,delta,-1],[0,-delta,-1]]
    aliceEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/4, createQubitObservable(bloch)),tetrahedronAlice))]
    
    
    bobUnBlochVectors = [[-1,-1,0],[-1,1,0]]
    bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    trineBob = [[0,0,-1],[-np.sqrt(3),0,1],[np.sqrt(3),0,1]]
    bobEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/3, createQubitObservable(bloch)),trineBob))]
    
    tetrahedronBob = [[1,-delta,0],[1,delta,0],[-1,0,delta],[-1,0,-delta]]
    tetrahedronBob = tetrahedronAlice
    bobEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/4, createQubitObservable(bloch)),tetrahedronBob))]
    
    aliceEffects = addEffectsForAbortOutcomes(aliceEffects,2)
    bobEffects = addEffectsForAbortOutcomes(bobEffects,2)
    
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
    
        
    
