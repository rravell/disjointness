import numpy as np
import cvxopt as cvx
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt


if __name__ == '__main__':

  outputsAlice = [4,4,4,4]
  outputsBob = [4,4,4,4]
  inputs=len(outputsAlice)
  vertices=generateLocalVertices(outputsAlice,outputsBob)
  symmetricVertices=Symmetrise(vertices,inputs,outputsAlice,outputsBob)
  print(symmetricVertices)
  
  
    
  
