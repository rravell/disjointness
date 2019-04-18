'''
Created on Apr 18, 2019

@author: rravell
'''


import numpy as np
N=4
K=2+1
monomialexpectedvalue=0.3

values={
    }

qdisjointness=np.zeros((N*K)**2)
for i in range (0, np.size(qdisjointness)):
    NumberOfGroup=i//(K*N)
    NumberOfCoefficient=i%(K*N)
    NumberOfInputA=NumberOfGroup//N 
    NumberOfInputB=NumberOfGroup%N                 
    NumberOfOutputA=NumberOfCoefficient//K 
    NumberOfOutputB=NumberOfCoefficient%K
    if NumberOfCoefficient<(K*(K-1)):
        if (i+1)%K==0:
            qdisjointness[i]=monomialexpectedvalue-qdisjointness[i-1]-qdisjointness[i-2]
        else:
            qdisjointness[i]=values['A{}|{} B{}|{}' .format(NumberOfOutputA, NumberOfInputA, NumberOfOutputB, NumberOfInputB)]
    else:
        if NumberOfCoefficient!=K**2-1:
            qdisjointness[i]=monomialexpectedvalue-qdisjointness[i-K] - qdisjointness[i-2*K]
        else:
            qdisjointness[i]=1-2*monomialexpectedvalue-qdisjointness[i-1]-qdisjointness[i-2]
            
print(qdisjointness)