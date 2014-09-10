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

def import_otu_header_files(file_1, file_2):
    karen_file = open(file_1, 'r')
    mark_file = open(file_2, 'r')

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
        
    karen_file.close()
    mark_file.close()
        
    return karen, mark

####





def main():
    if len(sys.argv) > 1: 
        two_files = sys.argv[1,2] # may need to separate into two arguments
    else:
        two_files = "karen.csv", "mark.csv"
        
    karen_names, mark_names = import_otu_header_files(two_files)

if __name__ == '__main__':
    main()