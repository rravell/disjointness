'''
Created on Apr 12, 2019

@author: rravell
'''


import numpy as np


N=4
K=2+1

#P y monomialexpectedvalue obtenidos en Lisp. 
p=0.413
monomialexpectedvalue=0.3287
#Distribucion de EqualityPrime clasica
vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #HAMMING DISTANCE N/2
vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #IGUALES
vector3=np.array([0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0]) #CASO CONTRARIO
EqualityPrime=np.concatenate((vector2,vector1,vector1,vector3,vector1,vector2,vector3,vector1,vector1,vector3,vector2,vector1,vector3,vector1,vector1,vector2))

#Obtencion de las componentes de qdisjoint para los eventos que no abortan
EqualityPrimeaborting=np.zeros((N*K)**2)
for i in range(0,np.size(EqualityPrimeaborting)):
    EqualityPrimeaborting[i]=EqualityPrime[i]*p


qEqualityPrime=np.array(EqualityPrimeaborting)


#Bucle para obtener qdisjoint para los eventos que si abortan
for i in range (0, np.size(qEqualityPrime)):
    NumberOfGroup=i//(K**2)
    NumberOfCoefficient=i%(K**2)
    if NumberOfCoefficient<(K*(K-1)):
        if (i+1)%K==0:
            qEqualityPrime[i]=monomialexpectedvalue-qEqualityPrime[i-1]-qEqualityPrime[i-2]
        else:
            continue
    else:
        if NumberOfCoefficient!=K**2-1:
            qEqualityPrime[i]=monomialexpectedvalue-qEqualityPrime[i-K] - qEqualityPrime[i-2*K]
        else:
            qEqualityPrime[i]=1-2*monomialexpectedvalue-qEqualityPrime[i-1]-qEqualityPrime[i-2]




print(qEqualityPrime)
np.savetxt('qEqualityPrimedist',qEqualityPrime)