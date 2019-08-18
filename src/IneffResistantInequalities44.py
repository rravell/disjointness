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
from bellscenario import BellScenario
from itertools import product
from bellpolytope import BellPolytope

if __name__ == '__main__':

    outputsAlice = [4,4,4]
    outputsBob = [4,4,4]
    K=outputsAlice[0]
    N=len(outputsAlice)
    inputs=len(outputsAlice)
    SymmetricVertices=generateLocalSymmetricVertices(outputsAlice, outputsBob, inputs)
    inequalities=np.loadtxt('Facets3-4.txt')
    functionals=-inequalities[:,1:]
    ineffFunctionals=[ToNormalConfiguration(functional, N, K) for functional in functionals]
    scenario=BellScenario(outputsAlice,outputsBob)
    poly=BellPolytope(scenario)
    vertices=[ver.getProbabilityList() for ver in poly.getListOfVertices()]
    normInef=list(map(lambda fun : 1/max([np.dot(fun,ver) for ver in vertices])*fun,ineffFunctionals))
    trueIneff = filter(lambda func :
                        not reduce(lambda coef,acum : coef==0 and acum,
                                   [func[scenario.getTuplesOfEvents().index((x,y,a,b))] 
                                    for x,y in product([0,1,2],repeat=2) for a,b in product([0,1,2,3],repeat=2) if a==2 or b==2 ],True),normInef)
    text=open("IneffFunctionals3-4Lisp.txt",'w')
    lisp=LispHelper
    for func in ineffFunctionals:
        lispsyntax=lisp.syntaxlisp(lisp, func, N, K)
        text.write(lispsyntax + '\n')
    np.savetxt('IneffFunctionals3-4.txt',ineffFunctionals)
    text.close()
