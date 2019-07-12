
'''
Created on Feb 12, 2019

@author: rravell
'''
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt

if __name__ == '__main__':
    
    outputsAlice = [3,3,3,3,3,3,5]
    outputsBob = [3,3,3,3,3,3,5]
    
    aliceUnBlochVectors = [[1,0,0],[0,1,0],[0,0,1]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceUnBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
    
    delta=0.01
    tetrahedronAlice = [[delta,0,1],[-delta,0,1],[0,delta,-1],[0,-delta,-1]]
    aliceEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/4, createQubitObservable(bloch)),tetrahedronAlice))]
    
    
    bobUnBlochVectors = [[1,0,0],[0,-1,0],[0,0,1]]
    bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    tetrahedronBob = [[1,-delta,0],[1,delta,0],[-1,0,delta],[-1,0,-delta]]
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
    
        
    
