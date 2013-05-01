# Checks update numbers for ties -- verify use of rounded-to-50 updates with
# with ranksum tests
#
# BLW 
# 8-30-09
# 7-11-11
#
# Throughout, manual checking is done instead of using sets so that information
# about WHICH seeds have tied updates may be reported.

import argparse

parser = argparse.ArgumentParser(description = "Report tied approximate updates within and across experiments.")

parser.add_argument("-f1", "--filegroup_1", metavar = "FILE", required = True, nargs='*', help = "_taskinfo.dat files from the first experiment (one is typical)")
parser.add_argument("-f2", "--filegroup_2", metavar = "FILE", required = True, nargs = '*', help = "_taskinfo.dat files from the second experiment (one is typical)")

args = parser.parse_args()
args = vars(args)  # turn namespace into dictionary

filegroup_1 = args["filegroup_1"]
filegroup_2 = args["filegroup_2"]

SEED = 0
FIRST_APPROX = 1
FIXED_APPROX = 2
FIRST_REAL = 3
FIXED_REAL = 4



# Collect all the numbers from an experiment:
def make_update_list(filegroup):
	update_list = []
	
	for filename in filegroup:
		fd = open(filename, "r")
		
		for line in fd:
			line = line.strip()
			
			if len(line) == 0 or line[0] == "#":
				continue
				
			line = line.split()
			
			# approximations
			first_EQU = int(line[1])
			first_EQU -= (first_EQU % 50)
			fixed_EQU = int(line[2])
			fixed_EQU -= (fixed_EQU % 50)
			
			update_list.append([int(line[0]), first_EQU, fixed_EQU, int(line[1]), int(line[2])])
		
		fd.close()
		
	return update_list
		
		
# We want to know:
# Are there ties in the approximate updates?
# If so, they can be ignored if they are also ties in the real updates,
# so we want to know that too.
#
# I'm not feeling particularly imaginative at the moment,
# so this will be straightforward n^2.
def report_ties(seeds_1, approx_1, real_1, seeds_2, approx_2, real_2):
	
	for i in range(0, 100):
		for j in range(0, 100):
			if approx_1[i] == approx_2[j] and not approx_1[i] < 0:
				print "Seed", seeds_1[i], "and seed", seeds_2[j], "both use update", str(approx_1[i]) + "."
				print "  (True seeds are", real_1[i], "and", str(real_2[j]) + ".)"
	
	

# Check against self
def report_self_ties(seeds, approx, real):
	
	for i in range(0, 100):
		for j in range(i + 1, 100):
			if approx[i] == approx[j] and not approx[i] < 0:
				print "Seed", seeds[i], "and seed", seeds[j], "both use update", str(approx[i]) + "."
				print "  (True seeds are", real[i], "and", str(real[j]) + ".)"
		
		
# MAIN PROGRAM BEGIN:		
		
update_list_1 = make_update_list(filegroup_1)
update_list_2 = make_update_list(filegroup_2)



# Are there ties within an experiment?  
# It probably doesn't matter but it might be nice to know.

print "Checking for ties within experiments: Group 1, first_EQU..."

report_self_ties([x[SEED] for x in update_list_1], [x[FIRST_APPROX] for x in update_list_1], [x[FIRST_REAL] for x in update_list_1])

print
print "Checking for ties within experiments: Group 1, fixed_EQU..."

report_self_ties([x[SEED] for x in update_list_1], [x[FIXED_APPROX] for x in update_list_1], [x[FIXED_REAL] for x in update_list_1])

print
print "Checking for ties within experiments: Group 2, first_EQU..."

report_self_ties([x[SEED] for x in update_list_2], [x[FIRST_APPROX] for x in update_list_2], [x[FIRST_REAL] for x in update_list_2])

print
print "Checking for ties within experiments: Group 2, fixed_EQU..."

report_self_ties([x[SEED] for x in update_list_2], [x[FIXED_APPROX] for x in update_list_2], [x[FIXED_REAL] for x in update_list_2])

# What we really want to know: are ties possibly affecting ranking?
# That is, check for ties between experiments.

print
print "Checking for ties between experiments: first_EQU..."

report_ties([x[SEED] for x in update_list_1], [x[FIRST_APPROX] for x in update_list_1], [x[FIRST_REAL] for x in update_list_1], [x[SEED] for x in update_list_2], [x[FIRST_APPROX] for x in update_list_2], [x[FIRST_REAL] for x in update_list_2])

print
print "Checking for ties between experiments: fixed_EQU..."

report_ties([x[SEED] for x in update_list_1], [x[FIXED_APPROX] for x in update_list_1], [x[FIXED_REAL] for x in update_list_1], [x[SEED] for x in update_list_2], [x[FIXED_APPROX] for x in update_list_2], [x[FIXED_REAL] for x in update_list_2])
