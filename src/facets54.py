import numpy as np
import cvxopt as cvx
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt
import cdd
from sympy.combinatorics.generators import symmetric


if __name__ == '__main__':

    outputsAlice = [4,4,4,4,4]
    outputsBob = [4,4,4,4,4]
    inputs=len(outputsAlice)
    SymmetricVertices=generateLocalSymmetricVertices(outputsAlice, outputsBob, inputs)
    
    VRepresentation=np.ones((len(SymmetricVertices), len(SymmetricVertices[0])+1))
    
    for i in range (0, len(SymmetricVertices)):
        VRepresentation[i][1:]=SymmetricVertices[i][:]
    
    mat = cdd.Matrix(VRepresentation, number_type='fraction')
    mat.rep_type = cdd.RepType.GENERATOR
    poly=cdd.Polyhedron(mat)
    print(poly.get_inequalities())
  
