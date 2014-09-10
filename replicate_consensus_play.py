import numpy as np

x = ([[0, 1, 9, 1, 0, 0],
  [1, 0, 1, 0,  1, 1],
  [1, 1,  0, 1,  1, 1],
  [1,  0,  0,  9,  0, 1],
  [ 0,  1,  1,  1,  1, 1]])

a=np.array(x) # this will be an array of genotypes where each row is a locus
#               and each column is a replicate

y = ([[7,7],[8,8],[11,11],[12,12],[14,14]])

b=np.array(y)   


"""avg = a.sum(axis=1)/float(a.shape[1])   # this is a list of the average score

#but now I want to change anything over 0.5 to '1', less than 0.5 as '0' and
# anything equal to 0.5 as '9'. This will be the consensus genotype and '9'
# means it is ambiguous
print a
print avg
n=[]

for f in avg:
    if (f > 0.3) and (f < 0.7): # can also use assert_almost_equals() in nose
        l = 9
        #count these
    elif f > 0.5:
        l=1
    elif f < 0.5:
        l=0
    else:
        l = f
    n.append(l)
print n
"""