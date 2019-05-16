
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
     
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables),np.float16)
     
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
    
    
    mat = cdd.Matrix(matrixHRepresentation, number_type='fraction')
    mat.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(mat)
    vertices = poly.get_generators()
    print(vertices)

    
