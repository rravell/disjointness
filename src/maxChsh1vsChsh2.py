
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
from itertools import product

import picos as pic

from math import sqrt

def CHSH(A,B):
    #chsh bell op give observables
    return pic.kron(A[0],B[0])+pic.kron(A[0],B[1])+pic.kron(A[1],B[0]) \
            -pic.kron(A[1],B[1])

if __name__ == '__main__':
    
     
    outputsAlice = [4,4,4,4]
    outputsBob = [2,2]
    
           
    alpha=2.7
    
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    
    prob=pic.Problem() 
    
    A={}
    for x1,x2 in product(range(2),range(2)):
        for a1,a2 in product(range(2),range(2)):
            A[x1,x2,a1,a2]=prob.add_variable('A_{0}{1}{2}{3}'.format(x1,x2,a1,a2),
                                             (2,2),'hermitian')
            prob.add_constraint(A[x1,x2,a1,a2]>>0)
    
    for x1,x2 in product(range(2),range(2)):
        prob.add_constraint(sum(A[x1,x2,a1,a2] for a1 in range(2) 
                                               for a2 in range(2))==np.eye(2))
    
    z0=np.array([[1,0],[0,0]])
    z1=np.array([[0,0],[0,1]])
    x0=np.array([[1,1],[1,1]])/2
    x1=np.array([[1,-1],[-1,1]])/2
    
    
    B={}
    B[0,0]=pic.new_param('B_00',z0)
    B[0,1]=pic.new_param('B_01',z1)
    B[1,0]=pic.new_param('B_10',x0)
    B[1,1]=pic.new_param('B_11',x1)
    
    rho=pic.new_param('rho',np.outer([1,0,0,1],[1,0,0,1])/2)
    
    
    A1=[1/2*sum(A[x1,x2,a1,a2]*(-1)**a1 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x2 in range(2))
                                                    for x1 in range(2)]
    
    A2=[1/2*sum(A[x1,x2,a1,a2]*(-1)**a2 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x1 in range(2))
                                                    for x2 in range(2)]
    
    B1=[sum(B[y,b]*(-1)**b for b in range(2)) for y in range(2)]
    
    
    prob.add_constraint(pic.trace(CHSH(A1,B1)*rho)==alpha)
    
    prob.set_objective('max',
                       pic.trace(CHSH(A2,B1)*rho))
    
    prob.solve()
    
    dist=list(np.zeros((64,1)))
    for x1,x2,y in product(range(2),repeat=3):
        for a1,a2,b in product(range(2),repeat=3):
            dist[(b+2*a1+4*a2)+8*(y+2*x2+4*x1)]=(
                                 pic.trace(pic.kron(A[x1,x2,a1,a2],B[y,b])*rho)
                                                                    ).get_value()
              
    


    
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
        
    


