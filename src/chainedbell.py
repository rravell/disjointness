
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
import cvxopt as cvx
from linopttools import *
import qutip as qt


def createEffects(eigenvectorsPerInput):
    return list(map(lambda eigenvectors :  
                   list(map(lambda eigvector : eigvector*eigvector.dag(),eigenvectors)),eigenvectorsPerInput))

if __name__ == '__main__':
    
    dim=3
    totalOutputs=dim+1
    N=4
    outputsAlice = [totalOutputs for i in range(0,N)]
    outputsBob = [totalOutputs for i in range(0,N)]
    
    
    eigVectAlice= [[1/np.sqrt(dim)*sum([np.exp(2*np.pi*1j/dim*q*(r-(k-1/2)/N))*qt.basis(dim, q) for q in range(0,dim)])
                    for r in range(0,dim)] 
                    for k in range(1,N+1)]
    
    eigVectBob= [[1/np.sqrt(dim)*sum([np.exp(-2*np.pi*1j/dim*q*(r-l/N))*qt.basis(dim, q) for q in range(0,dim)])
                    for r in range(0,dim)] 
                    for l in range(1,N+1)]
    
    aliceEffects = addEffectsForAbortOutcomes(createEffects(eigVectAlice),dim)
    bobEffects = addEffectsForAbortOutcomes(createEffects(eigVectBob),dim)
    psi=createMaxEntState(dim)
     
    dist = computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
    dist = [p.real for p in dist]
     
    vertices = generateLocalVertices(outputsAlice, outputsBob)
    with Model("lo1") as M:

        # Create variable 'x' of length 4
        bellFunctional = M.variable("bellFunctional", len(vertices[0]))
        
        i=0
        for vertice in vertices:
            M.constraint('c'+str(i),Expr.dot(vertice,bellFunctional), Domain.lessThan(1))
            i+=1
        
        for x,y in it.product(range(N),repeat=2):
            for a,b in [(i,dim) for i in range(0,dim)]+[(dim,i) for i in range(0,dim)]+[(dim,dim)]:
                indexPointer=list(np.zeros(len(vertices[0])))
                index=(b+totalOutputs*a)+((totalOutputs)**2)*(y+N*x)
                indexPointer[index]=1
                M.constraint('p('+str(index)+')',Expr.dot(bellFunctional,indexPointer),Domain.equalsTo(dist[index]))
         
        
        
        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, Expr.dot(dist, bellFunctional))
    
        # Solve the problem
        M.solve()
    
        # Get the solution values
        print(M.getPrimalSolutionStatus())
        print(M.primalObjValue())
