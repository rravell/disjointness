
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
import cvxopt as cvx
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt
    
if __name__ == '__main__':
    
    outputsAlice = [3,3,3,3]
    outputsBob = [3,3,3,5]
    
    aliceBlochVectors = [[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
    
    bobObservables=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    bobEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/4, createQubitObservable(bloch)),aliceBlochVectors))]
    
    aliceEffects = addEffectsForAbortOutcomes(aliceEffects,2)
    bobEffects = addEffectsForAbortOutcomes(bobEffects,2)
    
    psi=createMaxEntState(2)
    
    vertices = generateLocalVertices(outputsAlice, outputsBob)
    dist=computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)

    with Model("lo1") as M:

        # Create variable 'x' of length 4
        bellFunctional = M.variable("bellFunctional", len(vertices[0]))
        
        i=0
        for vertice in vertices:
            M.constraint('c'+str(i),Expr.dot(vertice,bellFunctional), Domain.lessThan(1))
            i+=1
        i=0
        for x in range (0,len(outputsAlice)):
            for y in range (0,len(outputsBob)):
                for a in range (0,outputsAlice[x]):
                    for b in range (0,outputsBob[y]):
                        indexPointer=list(np.zeros(len(vertices[0])))
                        indexPointer[i]=1
                        if(a==outputsAlice[x]-1)&(b==outputsBob[y]-1):
                            M.constraint('p('+str(i)+')',Expr.dot(bellFunctional,indexPointer),Domain.equalsTo(dist[i]))
                        else:
                            if(a==outputsAlice[x]-1)|(b==outputsBob[y]-1):
                                M.constraint('p('+str(i)+')',Expr.dot(bellFunctional,indexPointer),Domain.equalsTo(dist[i]))
                        i+=1  
        
        
        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, Expr.dot(dist, bellFunctional))
    
        # Solve the problem
        M.solve()
    
        # Get the solution values
        print(M.getPrimalSolutionStatus())
        print(M.primalObjValue())
        
    
    
        
    
