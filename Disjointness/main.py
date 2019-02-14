'''
Created on Feb 12, 2019

@author: rravell
'''
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np


if __name__ == '__main__':
    inputs = 4
    outputs = 2
     
    K = outputs+1
    totalOutputs = K**2
     
    N=inputs
    totalInputs = N**2
     
    nZeroConstraints = 2*totalInputs*(totalOutputs-(K-1)**2)
    nLowerBoundOnDeterministicConstraints = (K**(2*N))
    nUpperBoundOnDeterministicConstraints = (K**(2*N))
    nConstraints =  nZeroConstraints + nLowerBoundOnDeterministicConstraints + nUpperBoundOnDeterministicConstraints
    nVariables = (K*N)**2
     
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables))
     
    #IMPOSING THAT THE BELL COEFFICIENTS FOR ABORT EVENTS ARE >= 0
    rowNum = 0
    for j in range (1, 1+K*(totalInputs)):
        matrixHRepresentation[j-1, K*j]=1
        rowNum+=1
    base=rowNum
    for j in range (0, totalInputs):
        matrixHRepresentation[base + j, (K-1)*K+1 + totalOutputs*j]= 1
        matrixHRepresentation[base + 1 + j, (K-1)*K+2 + totalOutputs*j] = 1
        rowNum+=2 
            
    
    poly = BellPolytope(N,K)
    #IMPOSING THAT THE BELL FUNCTIONAL IS LOWER BOUNDED BY N**2*(N**2-2) ON LOCAL DISTRIBUTIONS WITH ABORT
    for vector in poly.getVertices():
        matrixHRepresentation[rowNum,0]=N**2*(N**2-2)
        matrixHRepresentation[rowNum,1:]=vector
        rowNum+=1
    
    #IMPOSING THAT THE BELL FUNCTIONAL IS UPPER BOUNDED BY 1 ON LOCAL DISTRIBUTIONS WITH ABORT    
    for vector in poly.getVertices():
        matrixHRepresentation[rowNum,0]=1
        matrixHRepresentation[rowNum,1:]=-vector
        rowNum+=1
    
    #IMPOSING THAT THE BELL COEFFICIENTS FOR ABORT EVENTS ARE <= 0
    base=rowNum
    for j in range (1,1+K*(totalInputs)):
        matrixHRepresentation[base+j,K*j]=-1
        rowNum+=1
    base=rowNum
    for j in range (0, totalInputs):
        matrixHRepresentation[base + j, (K-1)*K+1 + totalOutputs*j]= -1
        matrixHRepresentation[base + 1 + j, (K-1)*K+2 + totalOutputs*j] = -1
        rowNum+=2                         
                                                   
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #NO DISYUNTOS
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #SI DISYUNTOS
    Disjoint=np.concatenate((vector2,vector2,vector2,vector2,vector2,vector2,vector1,vector1,vector2,vector1,vector2,vector1,vector2,vector1,vector1,vector1))
    
    
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.obj_type = cdd.LPObjType.MAX
    mat.obj_func = tuple(np.concatenate(([0],Disjoint)))
    
    lp = cdd.LinProg(mat)
    lp.solve()
    print(lp.status)
    print(lp.obj_value)
    print(" ".join("{0}".format(val) for val in lp.primal_solution))
    
        
    