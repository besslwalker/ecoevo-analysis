# Designed to be called from calc_generations.sh -- if not so called, 
# there is no guarantee that the necessary time_tails.dat will exist.
#
# Reads time_tails.dat and parses out the number of generations at the end 
# of each run.  Calculates and reports the 5th percentile, median, and 95th
# percentile number of generations.
#
# BLW
# 9-27-11

import scipy.stats as ss

fd = open("time_tails.dat", 'r')

generations = []

for line in fd:
	line = line.strip()
	line = line.split()
	generations.append(float(line[2]))
	
print "5th percentile: " + str(ss.scoreatpercentile(generations, 5))
print "50th percentile: " + str(ss.scoreatpercentile(generations, 50))
print "95th percentile: " + str(ss.scoreatpercentile(generations, 95))