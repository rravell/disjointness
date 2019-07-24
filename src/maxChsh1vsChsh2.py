
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
    
           
    alpha=2.4
    
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
                                                                    ).get_value().real
              
    


    
    vertices=generateVertices1bitOfCommLocalPol(outputsBob,outputsAlice) 
     
    with Model("lo1") as M:

        # Create variable 'x' of length 4
        x = M.variable("x", len(vertices[0])+1)

        # Create constraints
        i=0
        for vertice in vertices:
            M.constraint('c'+str(i),Expr.dot(vertice+[-1],x), Domain.lessThan(0))
            i+=1
        M.constraint('c'+str(i),Expr.dot(dist+[-1],x),Domain.lessThan(10))
        
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
        print(dist)
#     
        chsh1 = lambda p : sum([(-1)**(x1*y)*1/2*(-1)**(a1+b)*p[(b+2*a1+4*a2)+8*(y+2*x2+4*x1)] for x2 in (0,1) for a1 in (0,1) for a2 in (0,1) for b in (0,1) for x1 in (0,1) for y in (0,1)])
        chsh2 = lambda p : sum([(-1)**(x2*y)*1/2*(-1)**(a2+b)*p[(b+2*a1+4*a2)+8*(y+2*x2+4*x1)] for x2 in (0,1) for a1 in (0,1) for a2 in (0,1) for b in (0,1) for x1 in (0,1) for y in (0,1)])
        dist=[0.42677669506557236, 0.0732233042012671, 4.3641591592724546e-10, -7.41575716416436e-10, -7.415757164248952e-10, 4.3641591596266833e-10, 0.0732233042012671, 0.42677669506557236, 0.42677669506558236, 0.07322330420125711, 4.3641591591640696e-10, -7.415757164055976e-10, -7.415757164211942e-10, 4.3641591595896737e-10, 0.07322330420125711, 0.42677669506558236, 0.24998754831860043, 1.2450977096787092e-05, 0.1267642439829991, 0.1232357555080078, 0.12323575550802345, 0.12676424398298997, 1.2450977096737804e-05, 0.24998754831862513, 0.1232357554117168, 0.12676424388398042, 1.2450977106337469e-05, 0.24998754851390056, 0.24998754851390723, 1.2450977106170935e-05, 0.12676424388398932, 0.12323575541173255, 0.24998754864831754, 1.2450977112536205e-05, 0.12676424381585538, 0.12323575534545382, 0.12323575534543428, 0.126764243815864, 1.2450977112617108e-05, 0.2499875486482894, 0.12676424405112593, 0.12323575557430413, 0.2499875481842188, 1.2450977090384951e-05, 1.2450977090572302e-05, 0.24998754818420768, 0.12323575557428447, 0.12676424405111755, 0.42677669506558025, 0.07322330420125918, 4.364159158812486e-10, -7.415757164439284e-10, -7.415757164566172e-10, 4.364159158712033e-10, 0.07322330420125917, 0.42677669506558036, 0.07322330420126502, 0.4267766950655744, -7.415757164317682e-10, 4.364159158690884e-10, 4.3641591586089363e-10, -7.415757164463076e-10, 0.4267766950655745, 0.07322330420126508]
    


