#This here will take a standard output which has not been formatted, and convert it into a more readable thing for the software

import sys

#Read in the 
input_file = open(sys.argv[1], 'r')

#Create the file to which we want to output:
output_file = open('formatted_data.txt', 'w+')


#Put up the first two lines needed
output_file.write('0\t0\t0\t50\t50\t50\t100\t100\t100\t500\t600\t500\t\n\n')

i = 0
for line in input_file:
    output_file.write('T%d\t' %i + line)
    i = i + 1