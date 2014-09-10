reg = open('reg1.csv', 'r')
l=[]
for line in reg:
    l.append(line.split(','))
print len(l[0])
