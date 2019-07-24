
'''
Created on Feb 12, 2019

@author: rravell
'''
from bellpolytope import BellPolytope
from mosek.fusion import *
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt

if __name__ == '__main__':
    
     
    outputsAlice = [2,2,2,2,3]
    outputsBob = [2,2]
        
    aliceUnBlochVectors = [[1,0,0],[0,0,1],[0,0,1]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceUnBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
     
    plus=1/np.sqrt(2)*(qt.basis(2, 0)+qt.basis(2, 1))
    minus=1/np.sqrt(2)*(qt.basis(2, 0)-qt.basis(2, 1))
 
    epsilon=0.001
    kraussPlus = np.cos(epsilon)*plus*plus.dag()+np.sin(epsilon)*minus*minus.dag()
    kraussMinus =  -np.cos(epsilon)*minus*minus.dag()+np.sin(epsilon)*plus*plus.dag()
    aliceEffects.append([kraussPlus*kraussPlus.dag(),kraussMinus*kraussMinus.dag()])
     
    symm3outDirect = [[0,0,1],[np.sin(2*np.pi/3),0,np.cos(2*np.pi/3)],[np.sin(4*np.pi/3),0,np.cos(4*np.pi/3)]]
    symm3outPovm = list(map(lambda bloch : 
                           effectForQubitPovm(1/3, createQubitObservable(bloch)),symm3outDirect))
    aliceEffects.append(symm3outPovm)
 
    mu = np.arctan(np.sin(2*epsilon))
    bobUnBlochVectors = [[np.sin(mu),0,np.cos(mu)],[-np.sin(mu),0,np.cos(mu)]]
    bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    psi=createMaxEntState(2)
    dist=computeDistributionFromStateAndEffects(psi,bobEffects,aliceEffects)
    
    vertices=generateVertices1bitOfCommLocalPol(outputsBob,outputsAlice) 
     
    with Model("lo1") as M:

        # Create variable 'x' of length 4
        x = M.variable("x", len(vertices[0])+1)

        # Create constraints
        i=0
        for vertice in vertices:
            M.constraint('c'+str(i),Expr.dot(vertice+[-1],x), Domain.lessThan(0))
            i+=1
        M.constraint('c'+str(i),Expr.dot(dist+[-1],x),Domain.lessThan(1))
        
        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, Expr.dot(dist+[-1], x))

        # Solve the problem
        M.solve()

        # Get the solution values
        print(M.getPrimalSolutionStatus())
        print(M.primalObjValue())
        print('local bound='+str(x.level()[len(vertices[0])]))
        bellFunctional=x.level()[0:len(vertices[0])]
        print('quantum value='+str(np.dot(bellFunctional,dist)))
        print(max([np.dot(bellFunctional,deterministicStg) for deterministicStg in vertices]))
#     


