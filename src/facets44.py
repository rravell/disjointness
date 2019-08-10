import numpy as np
import cvxopt as cvx
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt


if __name__ == '__main__':

  outputsAlice = [2,2]
  outputsBob = [2,2]
  inputs=len(outputsAlice)
  vertices=generateLocalVertices(outputsAlice, outputsBob, inputs)
  print(vertices)
  
    
  
