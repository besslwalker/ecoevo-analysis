# Takes the information produced by lineage_test.sh from lineage_test_raw.dat 
# and manipulates it into a human-readable form so the human can easily check 
# any runs that have final EQU-performing genotypes that did not correctly 
# get a lineage calculated.
#
# It produces the human-readable file lineage_test_suspicious.dat
#
# The human should examine all the suspicious lineages by hand.
#
# BLW
# 6-27-11

try:
	fd = open("lineage_test_raw.dat", "r");
except IOError:
	print "Could not open lineage_test_raw.dat"
	system.exit()


# File is formatted with a pair of lines per run:
# num_lines file_name
# lineage_depth blah blah blah blah blah blah task.0 task.1 task.2 task.3 task.4 task.5 task.6 task.7 task.8

first_line = True # housekeeper, flags whether this is the first or second line

suspicious_list = []

for line in fd:
	line = line.strip()
	line = line.split()

	if first_line:
		# extract number of lines and file name
		num_lines = int(line[0])
		file_name = line[1]
		
		# process file_name to get seed
		# file_name is <run-name>_<seed>/data/lineage_domEQU.dat
		file_name = file_name.split("/")
		# keep just the <run-name>_<seed> portion
		run_name = file_name[0]
		run_name = run_name.split("_")
		# extract just the <seed> portion and int-ify it
		seed = int(run_name[1])
		
	else:
		# extract lineage depth and EQU-performance
		depth = int(line[0])
		does_equ = bool(int(line[16])) 
	
		# Now that we have all the information, is this a suspicious run?		
		
		# We only care if it was an EQU-performing run
		if does_equ:
			# There are 23 blank/comment lines in the lineage_domEQU.dat files
			# Does the depth correspond to the number of lines?
			
			if num_lines - depth != 23:
				run = {}
				run["SEED"] = seed
				run["LINES"] = num_lines
				run["DEPTH"] = depth
				suspicious_list.append(run)
	
	first_line = not first_line # switch housekeeper value
	
fd.close()

# Now output those suspicious runs
try:
	fd = open("lineage_test_suspicious.dat", "w")
except IOError:
	print "Could not open lienage_test_suspicious.dat"
	system.exit()
	
if suspicious_list:
	for run in suspicious_list:
		fd.write(str(run["SEED"]) + " " + str(run["LINES"]) + " " + str(run["DEPTH"]) + "\n")
	print "Suspicious run seeds have been saved in lineage_test_suspicious.dat"
else:
	fd.write("No runs were suspicious!")
	print "No suspicious runs!"

	
fd.close()