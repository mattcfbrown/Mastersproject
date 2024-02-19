file = open("/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Data/Test2.txt", 'r')

f = open("/Users/mbrown/Desktop/Research/Mastersproject/Pipelines/Data/Test2_new.txt", 'w')

#Heading
f.write('0\t50\t100\t200\t300\n\n')
lines = file.readlines()

for i in range( len(lines) ):
    f.write('G' + repr(i) + '\t' + lines[i] + '\n')