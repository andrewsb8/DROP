import sys
import math

rama_file = open(sys.argv[1], 'r')
rama_data = rama_file.readlines()
rama_file.close()

output_file = open(sys.argv[2], 'w')

kT = 2.479 #kJ/mol @ 298 K

#get minimum nonzero value
min = 10 #probabilities will be less than 1
for line in rama_data:
    line = line.strip().split()
    if line == []:
        continue
    if float(line[2]) < min and float(line[2]) > 0:
        min = float(line[2])

for line in rama_data:
    line2 = line.strip().split()
    if line2 == []:
        output_file.write("\n")
    elif float(line2[2]) == 0:
        output_file.write(line)
    else:
        energy = -kT*math.log(float(line2[2])/min)
        output_file.write(line2[0] + "\t" + line2[1] + "\t" + str(energy) + "\n")

output_file.close()
