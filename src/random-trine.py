
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

def createRandomTrine():
    phase=np.random.uniform(-np.pi/2,np.pi/2)
    vectors=[[np.sin(phase+k*2*np.pi/3),0,np.cos(phase+k*2*np.pi/3)] for k in range(0,3)]
    return list(map(lambda vector: effectForQubitPovm(1/3, createQubitObservable(vector)),vectors))

if __name__ == '__main__':
    
    outputsAlice = [4,4,4,4]
    outputsBob = [4,4,4,4]
    
    aliceTrines = [createRandomTrine() for _ in range(0,len(outputsAlice))]
    bobTrines = [createRandomTrine() for _ in range(0,len(outputsBob))]
        
    aliceEffects = aliceTrines
#     aliceEffects[0] = projectorsForQubitObservable(createQubitObservable([1,0,0]))
    aliceEffects = addEffectsForAbortOutcomes(aliceEffects,2)
#     
    bobEffects = bobTrines
#     bobEffects[0] = projectorsForQubitObservable(createQubitObservable([0,1,0]))
#     bobEffects[1] = projectorsForQubitObservable(createQubitObservable([0,0,1]))
    bobEffects = addEffectsForAbortOutcomes(bobEffects,2)
    psi=createMaxEntState(2)
    
    dist=computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)

    vertices = generateLocalVertices(outputsAlice, outputsBob)
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
        
    
