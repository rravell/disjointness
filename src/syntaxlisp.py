'''
Created on 27 dic. 2018

@author: rravell
'''
import numpy as np

class LispHelper:
    def __init__(self):
        pass
    
    def genEquation(self, functional, j, N, K):
        NumberOfVector=(j)//(K**2) 
        NumberOfCoefficient=(j)%(K**2) 
        NumberOfInputA=NumberOfVector//N 
        NumberOfInputB=NumberOfVector%N                 
        NumberOfOutputA=NumberOfCoefficient//K 
        NumberOfOutputB=NumberOfCoefficient%K
            
        Number=abs(functional[j])
        syntax= "{} A{}/{} B{}/{}" .format(Number, NumberOfOutputA, NumberOfInputA, NumberOfOutputB, NumberOfInputB)
          
        return syntax

        NumberOfInputB=NumberOfVector%N           
    
    def syntaxlisp(self,functional, N, K):
        j = 0
        if functional[j]>=0:
            lispCommand = self.genEquation(self, functional, j, N, K)
        else:
            lispCommand = " - " + self.genEquation(self, functional, j, N, K)
        for j in range(1, np.size(functional)):
            if functional[j]>=0:
                lispCommand += " + " + self.genEquation(self, functional, j, N, K)
            else:
                lispCommand += " - " + self.genEquation(self, functional, j, N, K)
        return lispCommand
                
                

