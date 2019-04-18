'''
Created on Apr 12, 2019

@author: rravell
'''

import numpy as np

data = open('qdisjoint2nditeration')
clean_data = []
for quantity in data:
    number_and_d0 = quantity.split(' ')[-1]
    only_number = float(number_and_d0[:-3])
    clean_data.append(only_number)


np.savetxt('expectedvaluesqdisjoint',clean_data)