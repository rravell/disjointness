'''
Created on 9 jul. 2019

@author: gsenno
'''
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
import qutip as qt

def createQubitObservable(unnormalizedBlochVector):
    normalizedBlochVector=1/np.linalg.norm(unnormalizedBlochVector)*np.array(unnormalizedBlochVector)
    paulies=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    return sum([paulies[i]*normalizedBlochVector[i] for i in range(0,3)])

def effectForQubitPovm(weight,qubitObservable):
    return weight*(qt.qeye(2)+qubitObservable)

def projectorsForQubitObservable(qubitObservable):
    return list(map(lambda sign : effectForQubitPovm(1/2, sign*qubitObservable),[1,-1]))


def addEffectsForAbortOutcomes(effects,dim):
    return list(map(lambda povm : povm+[0*qt.qeye(dim)],effects))

def createMaxEntState(dimension):
    return 1/np.sqrt(dimension)*(sum([qt.tensor(qt.basis(dimension,i),
                                                   qt.basis(dimension,i)) for i in range(0,dimension)]))

def computeDistributionFromStateAndEffects(state,aliceEffects,bobEffects):
    dist=[]
    for x in range (0,len(aliceEffects)):
        for y in range (0,len(bobEffects)):
            for a in range (0,len(aliceEffects[x])):
                for b in range (0,len(bobEffects[y])):
                    dist.append((qt.tensor(aliceEffects[x][a],bobEffects[y][b])*state*state.dag()).tr())
    return dist                

def createConstraintsMatrixForEfficiencyDual(outputsAlice, outputsBob):
    KAlice = np.sum(outputsAlice)
    KBob = np.sum(outputsBob)
    
    
    totalOutputs = np.sum([a*b for a in outputsAlice for b in outputsBob])
    abortOutputs = totalOutputs - np.sum([(a-1)*(b-1) for a in outputsAlice for b in outputsBob])
    
    NAlice=len(outputsAlice)
    NBob=len(outputsBob)
    totalInputs = NAlice*NBob
     
    nZeroConstraints = 2*abortOutputs
    #nLowerBoundOnDeterministicConstraints = np.prod(outputsAlice)*np.prod(outputsBob)
    nLowerBoundOnDeterministicConstraints = 0
    nUpperBoundOnDeterministicConstraints = np.prod(outputsAlice)*np.prod(outputsBob)
    nConstraints =  nZeroConstraints + nLowerBoundOnDeterministicConstraints + nUpperBoundOnDeterministicConstraints
    nVariables = totalOutputs
     
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables),np.float16)
     
    #IMPOSING THAT THE BELL COEFFICIENTS FOR ABORT EVENTS ARE == 0
    #BOB'S ABORTS
    rowNum=0
    for i in range (0,NAlice):
        for j in range (0,NBob):
            for k in range (0,outputsAlice[i]):
                bellcoeff={}
                bellcoeff[toString(k, outputsBob[j]-1, i, j)]=1
                matrixHRepresentation[rowNum,1:]=toVector(bellcoeff,outputsAlice,outputsBob)
                bellcoeff[toString(k,outputsBob[j]-1,i,j)]=-1
                matrixHRepresentation[rowNum+1,1:]=toVector(bellcoeff,outputsAlice,outputsBob)
                rowNum+=2
            for k in range (0,outputsBob[j]-1):
                bellcoeff={}
                bellcoeff[toString(outputsAlice[i]-1,k,i,j)]=1
                matrixHRepresentation[rowNum,1:]=toVector(bellcoeff,outputsAlice,outputsBob)
                bellcoeff[toString(outputsAlice[i]-1,k,i,j)]=-1
                matrixHRepresentation[rowNum+1,1:]=toVector(bellcoeff,outputsAlice,outputsBob)
                rowNum+=2
    
       
    #IMPOSING THAT THE BELL FUNCTIONAL IS LOWER BOUNDED BY N**2*(N**2-2) ON LOCAL DISTRIBUTIONS WITH ABORT
#     for vector in poly.getVertices():
#         matrixHRepresentation[rowNum,0]=N**2*(N**2-2)
#         matrixHRepresentation[rowNum,1:]=vector
#         rowNum+=1
     
    #IMPOSING THAT THE BELL FUNCTIONAL IS UPPER BOUNDED BY 1 ON LOCAL DISTRIBUTIONS WITH ABORT    
    aliceStrategies = generateStrategies(outputsAlice)
    bobStrategies = generateStrategies(outputsBob)
    
    vertices=[strategyToDistribution(stgAlice,stgBob,outputsAlice,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]
 
    #vertices=BellPolytope(4,3).getVertices()
    for vector in vertices:
        matrixHRepresentation[rowNum,0]=1
        matrixHRepresentation[rowNum,1:]=-np.array(vector)
        #matrixHRepresentation[rowNum+1,0]=4**2*(4**2-2)
        #matrixHRepresentation[rowNum+1,1:]=np.array(vector)
        rowNum+=1
    
    return matrixHRepresentation    
        

def toString(a,b,x,y):
    return str(a)+','+str(b)+'|'+str(x)+','+str(y)
    
def toVector(dictionary,outputsAlice,outputsBob):
    vector = []
    for x in range (0,len(outputsAlice)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAlice[x]):
                for b in range (0,outputsBob[y]):
                    if toString(a, b, x, y) in dictionary:
                        vector.append(dictionary[toString(a, b, x, y)])
                    else:
                        vector.append(0)
    return vector
    

def generateStrategies(outputs):
    l=list(it.product(list(range(0,max(outputs))),repeat=len(outputs)))
    return list(it.filterfalse(lambda tup : not reduce(lambda x,y : x&y,list(imap(lambda a,b : a<b,tup,outputs))),l))
    
def strategyToDistribution(stgAlice, stgBob,outputsAlice,outputsBob):
    vector = []
    for x in range (0,len(outputsAlice)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAlice[x]):
                for b in range (0,outputsBob[y]):
                    if (a==stgAlice[x])&(b==stgBob[y]):
                        vector.append(1)
                    else:
                        vector.append(0)
    return vector
