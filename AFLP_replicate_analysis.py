import csv
import numpy as np
import sys
"""
1. Curate_aflp checks to make sure no weird or missing data.
2. The use separate script (compare_2_files.py) to report differences
   between the two scorers.
3. Differences between scorers are checked manually and rescored into a 
   reconciled csv file.
4. Manually removed H2O lanes.
5. Replace “unknown” samples (blind replicates) with correct sample names using
   standardize_sample_names and make_blind_rep_key.
6. Check for a bad lanes (find_bad_lanes). Pull all replicates of same sample. 
   For each replicate score the match to others of the same sample. If 7 
   replicates (for example) with 6 zeros and one lane with a "1", the lane with
   "1" should get a high mismatch score. Sum this across loci and it should give
   an indication of a bad lane (it will have a high combined mismatch score.)
7. Then discard bad lanes (manually in reconciled cvs file)
8. Then take each reconciled combo, pulling out only samples that were repeated,
   calculate a reliability score based on how consistent a locus is across
   replicates (find_bad_loci). Then we can discard bad loci, manually in csv file.
9. Now need to calculate the actual error rates across replicates. Done by
   find_bad_loci.
10. Next we need to generate a single column (lane) for each sample, removing
    all the extra replicates. Use simple majority rule. If equal (ratio between
    0.3 and 0.7) code as '9' and count. Done using function
    make_file_replicates_reconciled (along with replace_reps_with_consensus)
11. Concatenate the spread sheets from each combo into a single matrix. This is
    done manually
12. Convert to GenAlex format using separate script convert_to_genAlex.py


"""



def curate_aflp(aflp_matrix):
    datafile = open(aflp_matrix)
    aflp_report_file = open("C_report.txt", 'w')
    reader = csv.reader(datafile)
    counter = 0
    sample_list = []
    aflp_table = []
    locus_list = []
    for row in reader:
        if counter ==0:
            sample_list = row # a list of all samples names
            del sample_list[0] #just removing empty cell from first column
            counter +=1
        else:
            col_numb = 0; new_row = []
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
                        report_row = "sample " + str(get_xl_col(col_numb+1)) + " empty for locus " + str(counter-1) + " (row " + str(counter) + "): " + "\n"
                        aflp_report_file.write(report_row)
                    else:
                        new_row.append(int('0'))# get error if not int!
                        report_row = "column " + str(get_xl_col(col_numb+1)) + " is nonbinary for locus " + str(counter-1) + " (row " + str(counter) + "): " + str(cell) + "\n"
                        aflp_report_file.write(report_row)
                    col_numb +=1
            aflp_table.append(new_row)
            report_row = "sample size for locus " + str(counter-1) + ": " + str(col_numb-1) + "\n"
            aflp_report_file.write(report_row)#important check of data entry
        counter +=1
    table_array = np.array(aflp_table)
    datafile.close()
    aflp_report_file.close()
    return aflp_report_file#, table_array, sample_list, locus_list

def get_xl_col(index):
    counter = 0
    a1 = " ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    a2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for first in a1:
        for second in a2:
            counter +=1
            if counter == (index):
                column = str(first) + str(second)
                break
    return column   

def get_consensus_rep_genotpe(rep_array):
    avg = rep_array.sum(axis=1)/float(rep_array.shape[1]) # this is a list of the average score
    consensus_column = []
    for locus_av in avg:
        if (locus_av > 0.3) and (locus_av < 0.7):
            consensus = 9 #change to '-1'?
        elif locus_av > 0.5:
            consensus = 1
        elif locus_av < 0.5:
            consensus = 0
        else:
            consensus = locus_av
            print "check data: non-binary entry"
        consensus_column.append(consensus)
    return consensus_column

def replace_reps_with_consensus(locus_list, replicate_list, table_array,
                                standard_sample_list):
    delete_index_list = []
    consensus_matrix = []
    rep_name_list = []
    for reps in replicate_list:
        new_reps = [x - 1 for x in reps] # subtract 1 from each index value
        #because table_array is indexed differently. Check everything later
        rep_array = table_array[:,new_reps]
        consensus_column = get_consensus_rep_genotpe(rep_array)
        consensus_matrix.append(consensus_column)
        rep_name_list.append(standard_sample_list[reps[0]])
        for rep in reps:
            delete_index_list.append(rep)
    consensus_array = np.array(consensus_matrix).T
    ambiguous_count = 0
    for line in consensus_matrix:
        for col in line:
            if col == 9:
                ambiguous_count+=1
    print '# of "9"s: ', ambiguous_count
    return consensus_array, rep_name_list, delete_index_list

def build_matrix_minus_reps(standard_sample_list, table_array, 
                            consensus_array, rep_name_list, delete_index_list):
    #first delete replicates columns from names and aflp data
    aflp_cols_to_delete = [x - 1 for x in delete_index_list] # for correct indexing of table_array
    aflp_minus_reps = np.delete(table_array, aflp_cols_to_delete, 1)
    #now convert names list to array to delete indices, then back to list
    names_array = np.array(standard_sample_list)
    del_names_array = np.delete(names_array, delete_index_list, 0)
    del_names_list = del_names_array.tolist()
    # Now add back new data
    new_names_list = del_names_list + rep_name_list
    new_aflp_array = np.hstack((aflp_minus_reps,consensus_array))
    return new_aflp_array, new_names_list
    
def make_file_replicates_reconciled(new_aflp_array, new_names_list, locus_list):
    minus_reps = open("C01_minus_replicates.csv", 'w')
    for name in new_names_list:
        minus_reps.write(name + ',')
    minus_reps.write('\n')
    #minus_reps.write('next row')
    for locus, data in zip(locus_list, new_aflp_array):
        minus_reps.write(locus + ',')
        for cell in data:
            minus_reps.write(str(cell) + ',')
        minus_reps.write('\n')
    minus_reps.close()
        



def count_blind_h2o(aflp_report_file, sample_list):
    blind_report = open("AFLP_replicate_report.txt", 'a')
    blind_count = 0
    for item in sample_list:
        if "ERIunk" in item:
            blind_count +=1
    row = "\n" + str(blind_count) + " blind samples"
    blind_report.write(row)
    h2o_count = 0
    for item in sample_list:
        if "H2O" in item:
            h2o_count +=1
    row = "\n" + str(h2o_count) + " water controls"
    blind_report.write(row)
    blind_report.close()



def find_replicate_indices(sample_list,sample):
    sample_array = np.array(sample_list)
    search_name = sample
    index = np.where(values == search_name)[0]
    return index

# to index array: table_array[:,0] returns first column

def make_index_of_replicate_names(sample_list):
    replicate_index_list = []
    for name in sample_list:
        index_list = []
        index = 0
        sample = name[0:5]
        for i in sample_list: 
            if sample in i:
                index_list.append(index)
            index +=1
        replicate_index_list.append(index_list)
# Now clean the list to get rid of repeated lists
    new_replicate_index_list = []
    all_numbers = []
    for item in replicate_index_list:
        if len(item) > 1:
            for number in item:
                if number not in all_numbers:
                    all_numbers.append(number)
                    if item not in new_replicate_index_list:
                        new_replicate_index_list.append(item)
    return new_replicate_index_list


def count_diffs_bet_reps(table_array, standard_samples_names, x, y): 
    #x=389; y=408 # now make Boolean list if identical values at each locus. 
    #Note that indices are correct for sample list but subract one in table array
    #print standard_samples_names[x], standard_samples_names[y]
    #print table_array[:,x-1]
    #print table_array[:,y-1]
    compare_samples = table_array[:,x-1] == table_array[:,y-1]
    count = 0
    for i in compare_samples:
        if not i:
            count +=1
    return count    #This function works well

def get_worst_lanes(replicate_list, sample_list, index):
    line = []
    for list in replicate_list:
        if index in list:
            line.append(sample_list[index])
            line.append(list)
    return line
    
def find_bad_lanes(locus_list, table_array, replicate_list, standard_sample_list):
    bad_lane_file = open("bad_lanes.csv", 'w')
    worst_lanes_report = open("worst_lanes.csv", 'w')
    worst_lane_list = []
    max_diff = 0
    bad_index = 0
    total_count = 1
    total_diffs = 0
    for list in replicate_list:
        ind_count = 1
        for index in list:
            if ind_count ==1:
                if max_diff >2:
                    line = get_worst_lanes(replicate_list, standard_sample_list, bad_index)
                    worst_lane_list.append(line)
                max_diff = 0
                bad_lane_file.write('\n' + str(standard_sample_list[index]) + ',' + str(index))
                first_index = int(index)
                ind_count+=1
            else: 
                diffs = count_diffs_bet_reps(table_array, standard_sample_list, first_index, index)
                bad_lane_file.write("," + str(diffs))
                total_count+=1
                total_diffs += diffs
                if diffs > max_diff:
                    max_diff = diffs
                    bad_index = index
                ind_count+=1
    #print total_count
    print "total = ", total_count * len(locus_list)
    print "total count: ", total_count
    #print len(replicate_list)
    print "# of errors = ", total_diffs
    print "error rate = ", (float(total_diffs)/(total_count * len(locus_list))) * 100, "%"
    for row in worst_lane_list:
        c=1
        name = ''
        line = ''
        for item in row:
            if c == 1:
                name = str(row[0])
                c+=1
            else:
                for col in item:
                    line+= name + ','
                    line+= str(get_xl_col(col+1)) + ','
                    line += str(table_array[:,col-1]) + ','
                    line+= '\n'
                c+=1
        line+= '\n'
        worst_lanes_report.write(line)
    bad_lane_file.close()
    worst_lanes_report.close()
    
def find_bad_loci(locus_list, table_array, replicate_list, standard_sample_list):
    locus_report = open("locus_report.txt", 'w')
    locus_index = 0
    total = 0
    for locus in locus_list:
        line = locus
        total_locus_diff = 0
        for list in replicate_list:
            c = 1
            for index in list:
                if c ==1:
                    samp1 = index
                    c+=1
                else:
                    samp2 = index
                    if samp1 >= 0 and samp2 >= 0 and samp1 < len(standard_sample_list) and samp2 < (len(standard_sample_list)-1):
                        locus_diff = abs(table_array[locus_index,samp1-1] - table_array[locus_index,samp2-1])
                        #print locus, locus_diff, samp1, samp2
                        c+=1
                        total_locus_diff+= locus_diff
            
        locus_index+=1
        total+= total_locus_diff
        #print locus, total_locus_diff
    #print ""
    #print len(standard_sample_list)
    #print "total =", total


def test_make_jaccard_matrix(table_array, sample_list): # this works to find columns
    c=1
    for column in table_array.T:
        if c<4:
            print column
        c+=1
        
        
def make_jaccard_matrix(table_array, sample_list):
    j_matrix = open("Jaccard_matrix.txt", 'w')
    j_matrix.write(','+ str(sample_list) + '\n') #first comma gives blank in first column for row headers
    j_list = []
    xcount = 1 #remove for full matrix
    for x_column in table_array.T: #iterate also along sample_list to get row headers
        zero_count = 0; diff = 0
        for y_column in table_array.T:
            for i,j in zip(x_column, y_column):
                if i == 0 and j == 0:
                    zero_count +=1
                else:
                    if i != j:
                        diff +=1
            jaccard = float(diff)/float(len(x_column) - zero_count)
            j_list.append(jaccard)
            zero_count = 0; diff = 0
        j_matrix.write(str(j_list) + "\n")
        j_list = []
        xcount +=1 
        if xcount%10 == 0: #  print-to-screen counter
            print xcount
    j_matrix.close()
    
def clean_up_matrix(textfile):# strip [ and quotes
    dirty_file = open(textfile, 'r')
    clean_matrix = open('Jaccard_matrix.csv', 'w')
    for line in dirty_file:
        newline = line.replace("'","").replace('[','').replace(']','')
        clean_matrix.write(newline)
    dirty_file.close()
    clean_matrix.close()

    
def compare_names_across_combos(file1,file2): #probably do this manually
    b=9
    # run standardize_sample_name_start first - just gets rid of hyphens
    #
    
def concatenate_AFLP_combos(table_arrays):
    x=2
    # might need another function to look for same names in differnet sheets
    #probably need to use pandas to do this


def make_blind_rep_key(blind_key_file):
    key = open(blind_key_file, 'r')
    blind_rep_key = []
    for line in key:
        newline = line.strip('\r\n').split(',')
        blind_rep_key.append(newline)
    #print blind_rep_key
    return blind_rep_key

def standardize_sample_names(path, blind_rep_key):
    datafile = open (path, 'r')
    reader = csv.reader(datafile)
    newfile = open('new.txt', 'w')
    aflp_table = []
    locus_list = []
    counter = 0
    sample_list = []
    #first read in sample names for header row
    for row in reader:
        if counter ==0:
            sample_list = row # a list of all samples names
        else:
            break
        counter+=1
    datafile.close()
    #close file then create standard sample list
    standard_sample_list = []
    new_name = ''
    for name in sample_list:
        if 'Locus' not in name and 'H20' not in name and 'H2O' not in name and 'ERIunk' not in name and 'ERIUNK' not in name:
            new_name = "s" + name[3:5] + name[6:8]
        elif 'ERIunk' in name or 'ERIUNK' in name:
            unk_name = "ERIunk_" + name[7:9]
            for line in blind_rep_key:
                if line[0] == unk_name[0:9]:
                    decoded = "s" + line[1]
            new_name = decoded
        else:
            new_name = name #This is just saving the name "locus" for first column
        standard_sample_list.append(new_name)
    #Now reopen file to read in data
    datafile = open (path, 'r')
    reader = csv.reader(datafile)
    newline = ''
    for name in standard_sample_list:
        newline += name + ','
    newline += '\n'
    newfile.write(newline)
    counter= 0
    for row in reader:
        aflp_row = []
        if counter >0:#ignore header
            col = 0
            for i in row:
                if col == 0:#read in locus name
                    locus_list.append(i)
                    newfile.write(i + ',')
                    col+=1
                else:
                    aflp_row.append(int(i))
                    newfile.write(i + ',')
                    col+=1
            aflp_table.append(aflp_row)
            newfile.write('\n')
        counter+=1
    newfile.close()
    datafile.close()
    table_array = np.array(aflp_table)
    #print table_array.shape
    return standard_sample_list, locus_list, table_array, sample_list


def main():
    if len(sys.argv) > 1: 
        path = sys.argv[1]
    else:
        path = "C09_minus_bad_loci.csv"
        blind_key_file = "mw_blind_key.csv"
    #   flp_report_file = curate_aflp(path) #turn on to test raw input from one scorer
    blind_rep_key = make_blind_rep_key(blind_key_file)
    standard_sample_list, locus_list, table_array, orig_samples  = standardize_sample_names(path, blind_rep_key)
    #count_blind_h2o(aflp_report_file, sample_list)
    replicate_list = make_index_of_replicate_names(standard_sample_list) #orig_samples for new files
    #find_bad_lanes(locus_list, table_array, replicate_list, standard_sample_list)#orig_samples for new files
    #find_bad_loci(locus_list, table_array, replicate_list, standard_sample_list)
    consensus_array, rep_name_list, delete_index_list = replace_reps_with_consensus(locus_list, replicate_list, table_array, standard_sample_list)
    new_aflp_array, new_names_list = build_matrix_minus_reps(standard_sample_list, table_array, consensus_array, rep_name_list, delete_index_list)
    make_file_replicates_reconciled(new_aflp_array, new_names_list, locus_list)
    #print locus_list
    #print replicate_list
    #print standard_sample_list
    #print count_diffs_bet_reps(table_array)
    #make_jaccard_matrix(table_array, sample_list)
    #clean_up_matrix('Jaccard_matrix.txt')
    #count_diffs_across_loci(table_array)
    #print table_array.shape
    
if __name__ == '__main__':
    main()