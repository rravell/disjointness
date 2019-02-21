
'''
Created on Feb 12, 2019

@author: rravell
'''
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np


if __name__ == '__main__':
    inputs = 2
    outputs = 2
     
    K = outputs+1
    totalOutputs = K**2
     
    N=inputs
    totalInputs = N**2
     
    nZeroConstraints = 2*totalInputs*(totalOutputs-(K-1)**2)
    nLowerBoundOnDeterministicConstraints = 0
    nUpperBoundOnDeterministicConstraints = (K**(2*N))
    nConstraints =  nZeroConstraints + nLowerBoundOnDeterministicConstraints + nUpperBoundOnDeterministicConstraints
    nVariables = (K*N)**2
     
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables))
     
    #IMPOSING THAT THE BELL COEFFICIENTS FOR ABORT EVENTS ARE == 0
    rowNum=0
    for j in range (1, K*(totalInputs)+1):
        matrixHRepresentation[rowNum, K*j]=1
        matrixHRepresentation[rowNum+1, K*j]=-1
        rowNum+=2
    
    for j in range (0, totalInputs):
        matrixHRepresentation[rowNum, (K-1)*K+1 + totalOutputs*j]= 1
        matrixHRepresentation[rowNum+1, (K-1)*K+1 + totalOutputs*j]= -1
        matrixHRepresentation[rowNum+2, (K-1)*K+2 + totalOutputs*j] = 1
        matrixHRepresentation[rowNum+3, (K-1)*K+2 + totalOutputs*j] = -1
        rowNum+=4 
            
    poly = BellPolytope(N,K)
#     #IMPOSING THAT THE BELL FUNCTIONAL IS LOWER BOUNDED BY N**2*(N**2-2) ON LOCAL DISTRIBUTIONS WITH ABORT
#     for vector in poly.getVertices():
#         matrixHRepresentation[rowNum,0]=N**2*(N**2-2)
#         matrixHRepresentation[rowNum,1:]=vector
#         rowNum+=1
    
    #IMPOSING THAT THE BELL FUNCTIONAL IS UPPER BOUNDED BY 1 ON LOCAL DISTRIBUTIONS WITH ABORT    
    for vector in poly.getVertices():
        matrixHRepresentation[rowNum,0]=1
        matrixHRepresentation[rowNum,1:]=-vector
        rowNum+=1
    
    #DISJOINTNESS                                              
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #NO DISYUNTOS
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #SI DISYUNTOS
    Disjoint=np.concatenate((vector2,vector,vector2,vector2,vector2,vector2,vector1,vector1,vector2,vector1,vector2,vector1,vector2,vector1,vector1,vector1))
    Disjoint2=np.concatenate((vector2,vector2,vector2,vector1))
    #EQ'
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #HAMMING DISTANCE N/2
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #IGUALES
    vector3=np.array([0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0]) #CASO CONTRARIO
    EqualityPrime=np.concatenate((vector2,vector1,vector1,vector3,vector1,vector2,vector3,vector1,vector1,vector3,vector2,vector1,vector3,vector1,vector1,vector2))
    
    #GHD
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #<X,Y> < -SQRT(N)
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #<X,Y> >= SQRT(N)
    vector3=np.array([0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0]) #CASO CONTRARIO
    GapHammingDistance=np.concatenate((vector2,vector3,vector3,vector1,vector3,vector2,vector1,vector3,vector3,vector1,vector2,vector3,vector1,vector3,vector3,vector2))
    
    
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = tuple(np.concatenate(([0],Disjoint2)))
    
    lp = cdd.LinProg(mat)
    lp.solve()
    print(lp.status)
    print(lp.obj_value)
    print(" ".join("{0}".format(val) for val in lp.primal_solution))
    
        
    