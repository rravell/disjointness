'''
Created on Feb 14, 2019

@author: ricard
'''

def lispsyntax(functional, N, K):
    text=open("lispsyntax.txt", 'wa')
    for j in range(0, np.size(functional)):
        if functional[j]!=0:
                NumberOfVector=(j)/(K**2) 
                NumberOfCoefficient=(j)%(K**2) 
                NumberOfInputA=NumberOfVector/N 
                NumberOfInputB=NumberOfVector%N                 
                NumberOfOutputA=NumberOfCoefficient/K 
                NumberOfOutputB=NumberOfCoefficient%K
            
                Number=functional[j]
                a= "{} A{}/{} B{}/{} + " .format(Number, NumberOfOutputA, NumberOfInputA, NumberOfOutputB, NumberOfInputB)
                text.write(a)
                
                

