import sys
import math

lj_file = open(sys.argv[1], 'r')
lj_data = lj_file.readlines()
lj_file.close()

peptide_water_file = open(sys.argv[2], 'r')
peptide_water_data = peptide_water_file.readlines()
peptide_water_file.close()

output_file = open(sys.argv[3], 'w')

for i in range(len(lj_data)):
    lj_line = lj_data[i].strip().split()
    pw_line = peptide_water_data[i].strip().split()
    if lj_line == []:
        output_file.write("\n")
    else:
        output_file.write(lj_line[0] + "\t" + lj_line[1] + "\t" + str(float(lj_line[2])+float(pw_line[2])) + "\n")

output_file.close()
