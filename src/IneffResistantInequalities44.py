'''
Created on Aug 11, 2019

@author: ricard
'''


import numpy as np
import itertools as it
from _functools import reduce
from linopttools import *
import qutip as qt
from sympy.combinatorics.generators import symmetric
from syntaxlisp import LispHelper

if __name__ == '__main__':

  outputsAlice = [4,4,4,4]
  outputsBob = [4,4,4,4]
  K=outputsAlice[0]
  N=len(outputsAlice)
  inputs=len(outputsAlice)
  SymmetricVertices=generateLocalSymmetricVertices(outputsAlice, outputsBob, inputs)
  inequalities=np.loadtxt('Facets44.txt')
  functionals=-inequalities[:,1:]
  size=np.shape(functionals)
  IneffFunctionals=np.zeros((size[0],(K*N)**2))
  text=open("IneffFunctionals44Lisp.txt",'w')
  lisp=LispHelper
  for i in range (0, size[0]):
      values=np.dot(SymmetricVertices,np.transpose(functionals[i]))
      maximum=np.amax(values)
      if maximum!=0:
          normalisation=np.array((1/maximum)*functionals[i])
          IneffFunctionals[i]=ToNormalConfiguration(normalisation, N, K)
      else:
          normalisation=np.array(functionals[i])
          IneffFunctionals[i]=ToNormalConfiguration(normalisation, N, K)
      lispsyntax=lisp.syntaxlisp(lisp, IneffFunctionals[i], N, K)
      text.write(lispsyntax + '\n')
  np.savetxt('IneffFunctionals.txt',IneffFunctionals)
  text.close()
