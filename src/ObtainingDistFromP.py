'''
Created on Apr 12, 2019

@author: rravell
'''


import numpy as np


N=4
K=2+1

#P y monomialexpectedvalue obtenidos en Lisp. 
p=0.2914
monomialexpectedvalue=0.2812
#Distribucion de disjointness clasica
vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #NO DISYUNTOS
vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #SI DISYUNTOS
Disjointness=np.concatenate((vector2,vector2,vector2,vector2,vector2,vector1,vector2,vector1,vector2,vector2,vector1,vector1,vector2,vector1,vector1,vector1))

#Obtencion de las componentes de qdisjoint para los eventos que no abortan
disjointnessaborting=np.zeros((N*K)**2)
for i in range(0,np.size(disjointnessaborting)):
    disjointnessaborting[i]=Disjointness[i]*p


qdisjointness=np.array(disjointnessaborting)


#Bucle para obtener qdisjoint para los eventos que si abortan
for i in range (0, np.size(qdisjointness)):
    NumberOfGroup=i//(K**2)
    NumberOfCoefficient=i%(K**2)
    if NumberOfCoefficient<(K*(K-1)):
        if (i+1)%K==0:
            qdisjointness[i]=monomialexpectedvalue-qdisjointness[i-1]-qdisjointness[i-2]
        else:
            continue
    else:
        if NumberOfCoefficient!=K**2-1:
            qdisjointness[i]=monomialexpectedvalue-qdisjointness[i-K] - qdisjointness[i-2*K]
        else:
            qdisjointness[i]=1-2*monomialexpectedvalue-qdisjointness[i-1]-qdisjointness[i-2]




print(qdisjointness)
np.savetxt('qdisjointnessdist',qdisjointness)