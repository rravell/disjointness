'''
Created on 9 jul. 2019

@author: gsenno
'''
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
import qutip as qt

def getFirstAliceMarginal(n,dist):
    pAlice1 = []
    for x1,y in it.product(range(n),repeat=2):
        for a1,b in it.product(range(2),repeat=2):
            pAlice1.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                for x2 in range(n)
                for a2 in range(2)]))
    return pAlice1

def getSecondAliceMarginal(n,dist):
    pAlice2 = []
    for x2,y in it.product(range(n),repeat=2):
        for a2,b in it.product(range(2),repeat=2):
            pAlice2.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                    for x1 in range(n)
                    for a1 in range(2)]))
    return pAlice2

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
 
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables),np.int8)
     
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
    #vertices=BellPolytope(4,3).getVertices()
    aliceStrategies = generateStrategies(outputsAlice)
    bobStrategies = generateStrategies(outputsBob)
    
    #vertices=[strategyToDistribution(stgAlice,stgBob,outputsAlice,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]
    for stgAlice in aliceStrategies:
        for stgBob in bobStrategies:
            vertex=strategyToDistribution(stgAlice,stgBob,outputsAlice,outputsBob)
            matrixHRepresentation[rowNum,0]=1
            matrixHRepresentation[rowNum,1:]=-np.array(vertex)
            rowNum+=1
            #matrixHRepresentation[rowNum+1,0]=4**2*(4**2-2)
            #matrixHRepresentation[rowNum+1,1:]=np.array(vector)
    
    return matrixHRepresentation    

def generateLocalSymmetricVertices(outputsAlice,outputsBob, inputs):
    aliceStrategies = generateStrategies(outputsAlice)
    bobStrategies = generateStrategies(outputsBob)
    
    return [strategyToDistributionSymmetric(stgAlice,stgBob, inputs, outputsAlice,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]
    
def generateLocalVertices(outputsAlice,outputsBob):    
    aliceStrategies = generateStrategies(outputsAlice)
    bobStrategies = generateStrategies(outputsBob)
    
    return [strategyToDistribution(stgAlice,stgBob,outputsAlice,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]

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
    return l
    #return list(it.filterfalse(lambda tup : not reduce(lambda x,y : x&y,list(imap(lambda a,b : a<b,tup,outputs))),l))
    
def strategyToDistributionSymmetric(stgAlice, stgBob, inputs, outputsAlice,outputsBob):
    vector = []
    for x in range (0,len(outputsAlice)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAlice[x]):
                for b in range (0,outputsBob[y]):
                    if (a==stgAlice[x])&(b==stgBob[y]):
                        vector.append(1)
                    else:
                        vector.append(0)
    return VerticesToCG(vector, inputs, outputsAlice, outputsBob)

def strategyToDistribution(stgAlice, stgBob, outputsAlice,outputsBob):
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

def VerticesToCG(vector, inputs, outputsAlice, outputsBob):
    vertice = []
    #Alice's marginals
    l=0
    for x in range (0,len(outputsAlice)):
        s=0
        for y in range (0,outputsBob[0]):
            s += vector[l+y]
        vertice.append(s)
        for a in range (0,len(outputsBob)):
            l+=outputsAlice[x]*outputsBob[a]
    
    #Bob's marginals
    l=0
    for z in range (0,len(outputsAlice)):
        s=0
        for t in range (0,outputsAlice[0]):
            s += vector[l+outputsBob[z]*t]
        vertice.append(s)
        l+=outputsAlice[0]*outputsBob[z]
        
    #The rest of probabilities
    s=0
    for w in range (0,len(outputsAlice)):
        for x in range (0,len(outputsBob)):
            vertice.append(vector[s])
            s+=outputsAlice[w]*outputsBob[x]
    return Symmetrise(vertice, inputs)

def Symmetrise(vector,inputs):
    SymmetricVertices = []
    #First I transform the vector to CG notation
    #vector=VerticesToCG(vertices,outputsAlice,outputsBob)
    #Then I fill the marginals
    for i in range (0,inputs):
        SymmetricVertices.append(1/2*(vector[i]+vector[i+inputs]))
    #I create the matrix with the rest of probabilities
    CoefficientMatrix=np.zeros((inputs,inputs))
    s=0
    for l in range (0,inputs):
        for w in range (0,inputs):
            CoefficientMatrix[l][w]=vector[s+2*inputs]
            s+=1
    #The rest of probabilities
    for i in range (0,inputs):
        SymmetricVertices.append(CoefficientMatrix[i][i])
    for i in range (0,inputs):
        for j in range (0,inputs):
            if (i<j):
                SymmetricVertices.append(1/2*(CoefficientMatrix[i][j]+CoefficientMatrix[j][i]))
    return SymmetricVertices


    
    
def generateVertices1bitOfCommLocalPol(outputsAlice,outputsBob):
    communicationStrgs=list(it.product([0,1],repeat=len(outputsAlice)))
    #strgsAlice = [[[stgAlice[i],comm[i]] for i in range(0,len(stgAlice))] 
    #              for stgAlice in generateStrategies(outputsAlice) for comm in communicationStrgs]
    strgsBob = generateStrategies(list(reduce(lambda acum,elem : acum+[elem,elem],outputsBob,[])))
    
    vertices = []
    for stgAlice in generateStrategies(outputsAlice):
        for comm in communicationStrgs:
            for stgBob in strgsBob:
                vector = []
                for x in range (0,len(outputsAlice)):
                    for y in range (0,len(outputsBob)):
                        for a in range (0,outputsAlice[x]):
                            for b in range (0,outputsBob[y]):
                                if (a==stgAlice[x])&(b==stgBob[2*y+comm[x]]):
                                    vector.append(1)
                                else:
                                    vector.append(0)
                vertices.append(vector)
    return vertices


def ToNormalConfiguration(functionals,inputs,outputs):
    N=inputs
    K=outputs
    newfunctionals=np.zeros((N*K)**2)
    for j in range (0, N):
        for l in range (0,K-1):#Tengo en cuenta el abort
            #Los de Alice
            newfunctionals[l+j*N*K**2]+=1/2*functionals[j]
            #Los de Bob
            newfunctionals[(l*K)+j*K**2]+=1/2*functionals[j]
    CoefficientMatrix=np.zeros((N,N))
    s=0       
    r=0
    for w in range (0, N):
        for j in range (0, N):
            if (w==j):
                CoefficientMatrix[w][j]=functionals[s+N]
                s+=1
            elif (w<j):
                CoefficientMatrix[w][j]=1/2*functionals[r+2*N]
                r+=1
    for w in range (0,N):
        for j in range (0,N):
            CoefficientMatrix[j][w]=CoefficientMatrix[w][j]
    
    counting=0
    for w in range (0,N):
        for j in range (0,N):
            newfunctionals[counting*K**2]+=CoefficientMatrix[j][w]
            counting+=1
        
     
        
    return newfunctionals
    


    

