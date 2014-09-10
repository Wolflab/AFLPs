import csv
import numpy as np

# First read in data from csv file, curate (checking for '0' and '1')
# convert '01' to '1' and 'null' to '0'
datafile = open("Combo8_test_18Apr13_scoresKM.csv")
reader = csv.reader(datafile)

counter = 0
sample_list = []
table = []
locus_list = []
for row in reader:
    if counter ==0:
        sample_list = row
        del sample_list[0] #just removing empty cell from first column
        counter +=1
    else:
        col_numb = 0
        new_row = []
        for cell in row:
            if col_numb ==0:
                locus_list.append(cell)
                col_numb =+1
            else:
                if cell == '0' or cell == '1':
                    new_row.append(int(cell))
                elif cell == "null":
                    new_row.append(int('0'))
                elif cell == "01":
                    new_row.append(int('1'))
                elif cell == "":
                    new_row.append(int('0'))
                    print "sample number", col_numb, "empty: ", cell                   
                else:
                    new_row.append(int(cell))
                    #print "sample number", col_numb, "has nonbinary data: ", cell
                col_numb +=1
            
        table.append(new_row)
        counter +=1
        print "sample size for locus", counter, ":", col_numb-1 #important check of data entry

#Now count blind and control samples
"""blind_count = 0
for item in sample_list:
    if "ERIunk" in item:
        blind_count +=1
print blind_count, "samples"
h2o_count = 0
for item in sample_list:
    if "H2O" in item:
        h2o_count +=1
print h2o_count, "samples"  """

table_array = np.array(table)# that works!
# to index array: table_array[:,0] returns first column
# before you do anything else, copy file and make into functions!

#test to find index of sample name
il = []
index = 0
for i in sample_list:
    if i == "ERI05_17R_EacgMaga_5Dec12":
        il.append(index)
    index +=1
print il

x=121; y=122
# to compare replicates:
comp_list = table_array[:,x] == table_array[:,y]
for i,locus in zip(comp_list,locus_list):
    if i == False:
        print sample_list[x], "and", y, "differ", "at", locus
        
# comp_list is list of Boolean comparisons between columns\, one per locus
        

#print sample_list
#print len(locus_list), "loci"
#convert list to numpy array
#print locus_list
#print table[0:2]

datafile.close()


        
        