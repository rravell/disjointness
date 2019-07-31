
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
import cvxopt as cvx
import itertools as it
from ncpol2sdpa.sdp_relaxation import imap
from linopttools import *
import qutip as qt
from itertools import product
import picos as pic



def CHAINED(n,A,B):
    result = 0
    for i in range(0,n-1):
        result+=pic.kron(A[i],B[i])+pic.kron(A[i+1],B[i])
    result+=pic.kron(A[n-1],B[n-1])-pic.kron(A[0],B[n-1])
    return result

def chainedBellValue(n,p):
    result = 0
    for i in range(n-1):
        for a,b in product(range(2),repeat=2):
            result+=(-1)**(a+b)*(p[(b+2*a)+4*(i+n*i)]+p[(b+2*a)+4*(i+n*(i+1))])
    for a,b in product(range(2),repeat=2):
        result+=(-1)**(a+b)*(p[(b+2*a)+4*(n-1+n*(n-1))]-p[(b+2*a)+4*(n-1+n*(0))])
    return result

def getFirstAliceMarginal(n,dist):
    pAlice1 = []
    for x1,y in product(range(n),repeat=2):
        for a1,b in product(range(2),repeat=2):
            pAlice1.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                for x2 in range(n)
                for a2 in range(2)]))
    return pAlice1

def getSecondAliceMarginal(n,dist):
    pAlice2 = []
    for x2,y in product(range(n),repeat=2):
        for a2,b in product(range(2),repeat=2):
            pAlice2.append(sum([1/n*dist[(b+2*a1+4*a2)+8*(y+n*x2+(n**2)*x1)]
                    for x1 in range(n)
                    for a1 in range(2)]))
    return pAlice2


if __name__ == '__main__':
    
    n=3
    outputsAlice = [4 for _ in range(0,n**2)]
    outputsBob = [2 for _ in range(0,n)]
    
    communicationStrgs=list(it.product([0,1],repeat=len(outputsAlice)))
    
    maxChainedValues=set([])
    
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
                pAlice1 = getFirstAliceMarginal(n, vector)
                pAlice2 = getSecondAliceMarginal(n, vector)
                print('('+str(chainedBellValue(n, pAlice1))+','+str(chainedBellValue(n, pAlice2))+')')
                maxChainedValues.add((chainedBellValue(n,pAlice1),chainedBellValue(n,pAlice2)))
    print(maxChainedValues)
    

