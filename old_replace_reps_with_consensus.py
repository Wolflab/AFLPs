

    
def old_replace_reps_with_consensus(path): #This is a more simple way of importing
    #data, but I probably will not use it/
    datafile = open (path, 'r')
    newfile = open('reps_removed.txt', 'w')
    aflp_table = []
    locus_list = []
    counter = 0
    sample_list = []
    for row in datafile:
        col_count = 0
        if counter == 0:
            sample_list = row
            counter+=1
        else:
            for col in row:
                if col_count ==0:
                    locus_list.append(col)
                    col_count+=1
                else:
                    aflp_table.append(int(col))
                    col_count+=1
            counter+=1
                    
            
            