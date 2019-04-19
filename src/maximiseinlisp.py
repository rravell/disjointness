'''
Created on Mar 7, 2019

@author: rravell
'''
import numpy as np


if __name__ == '__main__':
    N=4
    K=2+1
    vector1=np.array([0.5, 0, 0, 0, 0.5, 0, 0 ,0 ,0]) #HAMMING DISTANCE N/2
    vector2=np.array([0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]) #IGUALES
    vector3=np.array([0.25, 0.25, 0, 0.25, 0.25, 0, 0, 0, 0]) #CASO CONTRARIO
    EqualityPrime=np.concatenate((vector2,vector1,vector1,vector3,vector1,vector2,vector3,vector1,vector1,vector3,vector2,vector1,vector3,vector1,vector1,vector2))
    text=open("maximiseinlisp.txt", 'w')
    for j in range (0, (N*K)**2):
        NumberOfVector=(j)//(K**2) 
        NumberOfCoefficient=(j)%(K**2) 
        NumberOfInputA=NumberOfVector//N 
        NumberOfInputB=NumberOfVector%N                 
        NumberOfOutputA=NumberOfCoefficient//K 
        NumberOfOutputB=NumberOfCoefficient%K
        Number=EqualityPrime[j]
        if NumberOfOutputB!=2 and NumberOfOutputA!=2:
            syntax= "(A{}/{} B{}/{} = {} p)" .format(NumberOfOutputA, NumberOfInputA, NumberOfOutputB, NumberOfInputB, Number)
            text.write(syntax +'\n' ) 
    text.close()