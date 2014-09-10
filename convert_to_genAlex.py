"""
this script takes output from AFLP_replcate_analysis.py. After manually
(in Excel) transposing data from each combo (selective primer set) and
concatanating all combos, this converts to genAlex format

"""
data_file = open("eri_aflp_all_combos.csv", "r")
genalex_file = open("eri_genalex.csv", 'w')
genalex_list = []
for line in data_file:
    genalex_list.append(line.strip().split(','))

#now count number of loci amd make list of ind names
number_of_inds = 0
ind_list = []
count = 1
for line in genalex_list:
    if count ==1: #count numnber of loci
        number_of_loci = len(line)-1
        count+=1
    else:
        number_of_inds+=1
        ind_list.append(line[0])# put pop names in list

#now use list of ind names to make list of pop names:            
pop_list = []
for ind in ind_list:
    pop = ind[0:3]
    if pop not in pop_list:
        pop_list.append(pop)

#now count inds in each pop
inds_in_pop_list = []  
count = 0
for pop in pop_list:
    inds_in_pop = 0
    for ind in ind_list:
        if pop in ind:
            inds_in_pop+=1
    inds_in_pop_list.append(inds_in_pop)
 
#now make line 1   
line1 = str(number_of_loci) + ',' + str(number_of_inds) + ',' + str(len(pop_list))
line = ''
for numb in inds_in_pop_list:
    line += ',' + str(numb)
line1+= line
line1+= '\n'

#now line 2:
line2 = 2*','
for name in pop_list:
    line2+= ',' + 'p' + name[1:]
line2+= '\n'
line3 = ','
    
genalex_file.write(line1)
genalex_file.write(line2)

#del genalex_list[0]#remove header before using
count = 1
for line in genalex_list:
    if count ==1:# get locus names
        newline = ','
        loc_count = 1
        for locus in line:
            if loc_count > 1: #to ignore the first item: 'locus' (the word)
                newline+= ',' + str(locus)
            else:
                loc_count+=1
        newline+= '\n'
        genalex_file.write(newline)
        count+=1
    else:
        sample_name = line[0]
        pop_name = 'p' + sample_name[1:3]
        newline = sample_name + ',' + pop_name
        
        col_count = 1
        for col in line:
            if col_count > 1:
                newline+= ',' + col
            else:
                col_count+=1
        newline+= '\n'
        genalex_file.write(newline)

      

data_file.close()
genalex_file.close()
    