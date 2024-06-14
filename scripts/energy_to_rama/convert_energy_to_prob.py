import sys
import math

rama_file = open(sys.argv[1], 'r')
rama_data = rama_file.readlines()
rama_file.close()

output_file = open(sys.argv[2], 'w')

kT = 2.479 #kJ/mol @ 298 K
sum = 0

#get min energy value
min = 1000000000000000 #probabilities will be less than 1
for line in rama_data:
    line = line.strip().split()
    if line == []:
        continue
    if float(line[2]) < min:
        min = float(line[2])

print(min)

for line in rama_data:
    line2 = line.strip().split()
    if line2 == []:
        continue
    else:
        try: #avoid overflow
            prob = math.exp(-( (float(line2[2]) - min) / kT))
            sum += prob
        except:
            sum += 0

print(sum)

#write normalized output
for line in rama_data:
    line2 = line.strip().split()
    if line2 == []:
        output_file.write("\n")
    else:
        try: #avoid overflow
            prob = math.exp(-( (float(line2[2]) - min) / kT))/sum
            output_file.write(line2[0] + "\t" + line2[1] + "\t" + str(prob) + "\n")
        except:
            output_file.write(line2[0] + "\t" + line2[1] + "\t" + "0\n")

output_file.close()
