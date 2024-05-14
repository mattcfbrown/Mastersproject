# #This here will take a standard output which has not been formatted, and convert it into a more readable thing for the software

import sys
import numpy as np

input_file = np.genfromtxt(sys.argv[1],
                           delimiter= '\t')

m = int(sys.argv[2])
iter = int(sys.argv[3])


output_file = open('formatted_data_' + str(iter) + '.txt', 'w+')
output_file.write('0\t0\t0\t50\t50\t50\t100\t100\t100\t500\t600\t500\t\n\n')


for i in range(len(input_file)):
    line = '\t'.join(map(str,input_file[i][m*(iter-1):(m-1)*iter]))
    output_file.write('T%d\t' %i + line)
    output_file.write('\n')