
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
import cvxopt as cvx
import itertools as it
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt
from itertools import product
import picos as pic
from gpg.results import Result

def CHAINED(n,A,B):
    result = 0
    for i in range(0,n-1):
        result+=pic.kron(A[i],B[i])+pic.kron(A[i+1],B[i])
    result+=pic.kron(A[n-1],B[n-1])-pic.kron(A[0],B[n-1])
    return result

def chainedBellValue(n,p):
    result = 0
    for i in range(n-1):
        for a,b in product(range(2),repeat=2):
            result+=(-1)**(a+b)*(p[(b+2*a)+4*(i+n*i)]+p[(b+2*a)+4*(i+n*(i+1))])
    for a,b in product(range(2),repeat=2):
        result+=(-1)**(a+b)*(p[(b+2*a)+4*(n-1+n*(n-1))]-p[(b+2*a)+4*(n-1+n*(0))])
    return result

if __name__ == '__main__':
    
    n=3
    outputsAlice = [4 for _ in range(0,n**2)]
    outputsBob = [2 for _ in range(0,n)]
    
#vertices=generateVertices1bitOfCommLocalPol(outputsAlice,outputsBob) 
           
    alpha=5.1
     
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
     
    prob=pic.Problem() 
     
    A={}
    for x1,x2 in product(range(n),repeat=2):
        for a1,a2 in product(range(2),repeat=2):
            A[x1,x2,a1,a2]=prob.add_variable('A_{0}{1}{2}{3}'.format(x1,x2,a1,a2),
                                             (2,2),'hermitian')
            prob.add_constraint(A[x1,x2,a1,a2]>>0)
     
    for x1,x2 in product(range(n),repeat=2):
        prob.add_constraint(sum(A[x1,x2,a1,a2] for a1 in range(2) 
                                               for a2 in range(2))==np.eye(2))
    
    phi = lambda i : (i-1)*np.pi/n
    phiprime= lambda i : (2*i-1)*np.pi/(2*n) 
    bobsProjectors=[projectorsForQubitObservable
                    (createQubitObservable([np.sin(phiprime(i)),0,np.cos(phiprime(i))]))
                     for i in range(1,n+1)]
    bobsProjectors=[list(map(lambda qutipProj : qutipProj.get_data().toarray(),obs)) for obs in bobsProjectors]
     
     
    B={}
    for i in range(0,n):
        for j in (0,1):
            B[i,j]=pic.new_param('B_'+str(i)+str(j),bobsProjectors[i][j])
     
    rho=pic.new_param('rho',np.outer([1,0,0,1],[1,0,0,1])/2)
     
     
    A1=[1/n*sum(A[x1,x2,a1,a2]*(-1)**a1 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x2 in range(n))
                                                    for x1 in range(n)]
     
    A2=[1/n*sum(A[x1,x2,a1,a2]*(-1)**a2 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x1 in range(2))
                                                    for x2 in range(n)]
     
    B1=[sum(B[y,b]*(-1)**b for b in range(2)) for y in range(n)]
     
     
    prob.add_constraint(pic.trace(CHAINED(n,A1,B1)*rho)==alpha)
     
    prob.set_objective('max',
                       pic.trace(CHAINED(n,A2,B1)*rho))
     
    prob.solve()
#     
    dist=list(np.zeros((n*(n**2)*8,1)))
    for x1,x2,y in product(range(n),repeat=3):
        for a1,a2,b in product(range(2),repeat=3):
            try:
                dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]=(
                                 pic.trace(pic.kron(A[x1,x2,a1,a2],B[y,b])*rho)
                                                                    ).get_value().real
            except:
                print('index'+str(x1))
               
    pAlice1 = []
    for x1,y in product(range(n),repeat=2):
        for a1,b in product(range(2),repeat=2):
            pAlice1.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                for x2 in range(n)
                for a2 in range(2)]))
    pAlice2 = []
    for x2,y in product(range(n),repeat=2):
        for a2,b in product(range(2),repeat=2):
            pAlice2.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                    for x1 in range(n)
                    for a1 in range(2)]))
    
    
    print(chainedBellValue(n,pAlice1),chainedBellValue(n,pAlice2))