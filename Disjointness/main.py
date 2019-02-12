'''
Created on Feb 12, 2019

@author: rravell
'''
from bellpolytope import BellPolytope
import picos as pic
import numpy as np


if __name__ == '__main__':
    N=4
    K=2
    poly = BellPolytope(N,K+1)
    vertices=np.array(poly.getVertices())
    NumberOfCoefficients=(N*(K+1))**2
    shape=np.shape(vertices)
    NumberOfVertices=shape[0]
    NumberOfAborts=N*((K+1)*N + K*N)
    #NumberOfRows=NumberOfVertices + 2*NumberOfAborts
    A=np.zeros((NumberOfVertices, NumberOfCoefficients))
    #A[NumberOfVertices:NumberOfRows,0]=0
    A[:,:]=vertices
    coefficients=np.ones(NumberOfVertices)
    #l=0
    #s=0
    '''for j in range (NumberOfVertices, NumberOfVertices + NumberOfCoefficients/(K+1) + 1):
        A[j, (K+1)*l]=1
        l=l+1
      
    for j in range (0, N*(K+1)):
        A[NumberOfVertices + NumberOfCoefficients/(K+1) +K*j, NumberOfCoefficients/(N*(K+1))-2 + NumberOfCoefficients/(N*(K+1))*j]= 1
        A[NumberOfVertices + NumberOfCoefficients/(K+1) + 1 + K*j, NumberOfCoefficients/(N*(K+1))-1 + NumberOfCoefficients/(N*(K+1))*j] = 1 
 
    for j in range (NumberOfVertices + NumberOfCoefficients/(K+1) + 1,NumberOfVertices + 2*NumberOfCoefficients/(K+1) + 1):
        A[j,(K+1)*s]=-1
        s=s+1

    for j in range (0, N*(K+1)):
        A[NumberOfVertices + 2*NumberOfCoefficients/(K+1) +K*j, NumberOfCoefficients/(N*(K+1))-2 + NumberOfCoefficients/(N*(K+1))*j]= -1
        A[NumberOfVertices + 2*NumberOfCoefficients/(K+1) + 1 + K*j, NumberOfCoefficients/(N*(K+1))-1 + NumberOfCoefficients/(N*(K+1))*j] = -1 '''
    
    
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #NO DISYUNTOS
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #SI DISYUNTOS
    Disjoint=np.concatenate((vector2,vector2,vector2,vector2,vector2,vector2,vector1,vector1,vector2,vector1,vector2,vector1,vector2,vector1,vector1,vector1))
    
    
    P = pic.Problem()
    A = pic.new_param('A', A)
    Disjoint = pic.new_param('Disjoint', Disjoint)
    C= pic.new_param('C', coefficients)
    B = P.add_variable('B',NumberOfCoefficients)
    P.add_constraint(A*B < C)
    for j in range (0, NumberOfCoefficients/(K+1) + 1):
        P.add_constraint(B[(K+1)*j-1]==0)
    for j in range (0, N*N): #ESTOS CONSTRAINT DEPENDEN DE N, SI CAMBIO N HAY QUE CAMBIARLOS
        P.add_constraint(B[(K+1)*K + (K+1)*(K+1)*j]==0)
        P.add_constraint(B[(K+1)*(K+1) - 2 + (K+1)*(K+1)*j]==0)
    objective = np.dot(Disjoint.T,B)
    P.set_objective('max', objective)
    P.solve(verbose=0,solver='cvxopt')
    B_opt = B.value
    print(B_opt)
        
    