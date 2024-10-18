#This here will take a standard output which has not been formatted, and convert it into a more readable thing for the software

import sys

#Read in the data
input_file = open(sys.argv[1], 'r')

#Create the file to which we want to output:
# num_genes = sys.argv[2]
# name = 'formatted_data_' + num_genes + '.txt'
name = 'formatted_data.txt'
output_file = open(name, 'w+')


#Put up the first two lines needed
output_file.write('0\t0\t0\t50\t50\t50\t100\t100\t100\t500\t600\t500\t\n\n')

i = 0
for line in input_file:
    output_file.write('T%d\t' %i + line)
    i = i + 1