import numpy as np
import cvxopt as cvx
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt
import cdd
from sympy.combinatorics.generators import symmetric


if __name__ == '__main__':

  outputsAlice = [4,4,4]
  outputsBob = [4,4,4]
  inputs=len(outputsAlice)
  SymmetricVertices=generateLocalSymmetricVertices(outputsAlice, outputsBob, inputs)
  FilteredSymmetricVertices=[list(t) for t in set(tuple(element) for element in SymmetricVertices)]
  VRepresentation=np.ones((len(FilteredSymmetricVertices), len(FilteredSymmetricVertices[0])+1))
  
  for i in range (0, len(FilteredSymmetricVertices)):
      VRepresentation[i][1:]=FilteredSymmetricVertices[i][:]

  mat = cdd.Matrix(VRepresentation, number_type='fraction')
  mat.rep_type = cdd.RepType.GENERATOR
  poly=cdd.Polyhedron(mat)
  print(poly.get_inequalities())
