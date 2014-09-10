import numpy as np
import scipy as sp

data = np.genfromtxt('Combo8_test_18Apr13_scoresKM.csv', dtype=str, delimiter=',')

#put the first row and first column in their own arrays and remove them from
#the main array
indivID = data[0, 1:]
locus = data[1:, 0]
data = sp.delete(data, 0, 0)
data = sp.delete(data, 0, 1)
print locus