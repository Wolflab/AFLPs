"""This script should take two input files and compare them for differences, writing the differences to a third file and printing the total number of differences to the screen

IMPORTAMT - be sure to save Excel file as Windows Comma Separated!!!"""

import csv

def get_xl_col(index):
    column = ''
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

#import the data
karen_file = open('C10_EaggMatc_ME.csv', 'r')
mark_file = open('C10_EaggMatc_JW.csv', 'r')

# read into csv.reader
karen_csv = csv.reader(karen_file)
mark_csv= csv.reader(mark_file)

#read into list of lists
karen =[]
for line in karen_csv:
    karen.append(line)
    
mark =[]
for line in mark_csv:
    mark.append(line)


#create output:
report_file = open("report_differences_C10_JW_vs_ME_June_24.txt", 'w')
count = 0

"""report_line = "locus" + "\t" + "sample" + "\t" + "column" + "\n" +"\n"
report_file.write(report_line)"""

sample_list_m =[]
sample_list_k =[]
locus_count = -1
for i, (m_row,k_row )in enumerate(zip(mark,karen)):
    locus_m = ''
    ind_count = 1
    sample_count = 0
    if count == 0:
        sample_list_m = m_row
        sample_list_k = k_row
        name_count = -1
        for sample_name, (m_name,k_name) in enumerate(zip(sample_list_m, sample_list_k)):
            if m_name != k_name:
                report_line = "{compare names " + str(m_name) + '}\n'
                report_file.write(report_line)                
            name_count +=1
        count +=1
        #print sample_list (first row is all the sample names)
    else:
        null_count = 0
        if count ==1:
            report_line = "\n" + "locus" + "\t" + "sample" + "\t" + "column" + "\n" + "\n"
            report_file.write(report_line)
        locus_m = m_row[0]
        locus_k = k_row[0]
        if float(locus_k) != float(locus_m):
            report_line = "{compare locus " + locus_k + '}  '
            report_file.write(report_line)
        ind_count = 0
        for cell, (m_cell, k_cell) in enumerate(zip(m_row, k_row)):
            if m_cell == '0' and k_cell == "null":
                null_count +=1
            elif k_cell == '0' and m_cell == "null":
                null_count +=1
            elif m_cell != k_cell:
                count += 1
                sample_count +=1
                report_line = str(locus_m) + "\t" + str(sample_list_m[ind_count]) + "\t" + str(get_xl_col(ind_count+1)) + "\n"
                report_file.write(report_line)
            ind_count +=1
        report_line = "locus " + str(locus_m) + " error rate= " + str(float(sample_count)/float(ind_count)) +'\n'
        report_file.write(report_line)
    locus_count+=1

report_line = "\n" + "number of differences: " + str(count-1)
report_file.write(report_line)
report_line = "\n" + "number of cells compared: " + str(locus_count*name_count)
report_file.write(report_line)
report_line = "\n" + "error rate = : " + str(float(count-1)/(float(locus_count)*name_count))
report_file.write(report_line)
report_line = "\n" + "null count = : " + str(null_count)
report_file.write(report_line)

#close files
report_file.close()
karen_file.close()
mark_file.close()
