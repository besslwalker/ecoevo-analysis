# Reads tasks.dat and lineage-domEQU.dat files, produces files with lines in this format:
# 
# seed update_of_first_EQU update_of_fixed_EQU task_change_at_first_EQU task_change_at_fixed_EQU fitness_ratio_at_fixed_EQU
#
# If EQU does not appear or does not fix, update numbers are -1; 
# task change is 0, fitness change is -1 (always -1 if not using lineage stats)
#
# These files can be used to generate stats about time of EQU acquisition,
# tradeoffs at EQU acquisition, phenotypic diversity at EQU acqisition (by
# pinpointing the updates at which to examine phenotypic diversity).
#
# BLW
# August 2009
# 9-15-09
# 9-20-09
# 9-6-11

import gzip
from optparse import OptionParser
import os.path

# Set up parser
usage = """usage: %prog [options] output_prefix run_1 [run_2 ...]

Run names should not include the directory's trailing slash.
Seed numbers will be gleaned from the directory names of the runs.
Output file will be saved as <output_prefix>_taskinfo.dat."""

parser = OptionParser(usage)
parser.add_option("-l", "--lineage", action = "store_true", 
                  dest = "use_lineage", default = True, 
                  help = "use lineage file to determine EQU fixation stats")
parser.add_option("-t", "--tasks", action = "store_false", dest = "use_lineage",
                  help = "use tasks file to determine EQU fixation stats")

(options, args) = parser.parse_args()

# Go!
tasks_info = []

for run_directory in args[1:]:
	seed = int(run_directory.split("_")[-1])  # obtain the seed
	
	record = [seed]
	first_EQU = -1
	fixed_EQU = -1
	first_change = 0
	fixed_change = 0
	fixed_fitness = -1
	original_fixed_fitness = -1
	
	# First we run through tasks.dat
	# If determining fixation stats from tasks.dat, this is all we need to do
	#
	# This WILL count very-late-appearing EQU -- appearing in something that does
	# not die or reproduce before 100,000 updates -- as a real EQU fixation,
	# whether or not that EQU would have gone away after 100,000 updates.
	
	# try to open tasks.dat; if it doesn't exist, try to open tasks.dat.gz
	task_filename = run_directory + "/data/tasks.dat"
	try:
		tfd = open(task_filename, 'r')
	except (IOError):
		tfd = gzip.open(task_filename + ".gz")
	
	tasklines = tfd.readlines()
	
	for line_num, line in enumerate(tasklines):
		line = line.strip()
		
		if len(line) == 0 or line[0] == '#':
			continue
			
		line = line.split()
		
		# if we haven't yet seen EQU, record the update number and task change
		if first_EQU == -1 and int(line[9]) > 0:  
			first_EQU = int(line[0]) 
			
			prev_line = tasklines[line_num - 1].strip().split()
			first_change = sum([bool(int(x)) for x in line[1:10]])                  \
									 - sum([bool(int(x)) for x in prev_line[1:10]])
			
		# Even if we're using lineage for fixation stats, 
		# we need to know which runs fixed
		# So we need an initial value of fixed_EQU, but fixed_change does not
		# have to be calculated if we're not doing tasks.dat-based stats.
		# (And we can't know the ratio, so leave it at "N".)
			
		# if EQU is not currently fixed, record when it appears
		if fixed_EQU == -1 and int(line[9]) > 0:  
			fixed_EQU = int(line[0])
			if not options.use_lineage:
				prev_line = tasklines[line_num - 1].strip().split()
				fixed_change = sum([bool(int(x)) for x in line[1:10]])                \
										 - sum([bool(int(x)) for x in prev_line[1:10]])
				
		# if we think EQU has fixed but it's now gone, return to unfixed state
		if not fixed_EQU == -1 and int(line[9]) <= 0:
			fixed_EQU = -1
			fixed_change = 0
			
	# close file
	tfd.close()
	
	
	# If determining fixation stats from a lineage file, do so
	# This method will refuse to recognize any final-organism-only EQU
	# appearance as a real fixation, as there is no child to check for
	# whether EQU really fixed or just popped up.
	if options.use_lineage:
		# Did this seed fix EQU?
		if not fixed_EQU == -1:
			lineage_filename = run_directory + "/data/lineage_domEQU.dat"
			if not os.path.exists(lineage_filename):
				print "There is no lineage file for seed " + str(seed) + ".  Either you have not generated the lineage files or EQU appears in this population only as a stochastically executed task.  This population will be recorded as NOT having fixed EQU."
				fixed_EQU = -1
				continue
			else:	
				lfd = open(lineage_filename, "r")
			
			lineagelines = lfd.readlines()
			
			# reset fixed_EQU for fixation-finding
			fixed_EQU = -1
			
			for line_num, line in enumerate(lineagelines):
				line = line.strip()
				
				if len(line) == 0 or line[0] == "#":
					continue
				
				line = line.split()
				
				# If EQU is not currently fixed, record when it appears
				# Unless this is the very last organism in the lineage
				if fixed_EQU == -1 and int(line[16]) > 0 and line_num != len(lineagelines) - 1:
					fixed_EQU = int(line[7])
					prev_line = lineagelines[line_num - 1].strip().split()
					fixed_change = sum([int(x) for x in line[8:16]])             \
											 - sum([int(x) for x in prev_line[8:16]])
					fixed_fitness = float(line[6])
											 
				# If we though EQU had fixed but it's now gone, return to unfixed state
				if not fixed_EQU == -1 and int(line[16]) <= 0:
					fixed_EQU = -1
					fixed_change = 0
					fixed_fitness = -1
					
			lfd.close()
			
			# Get fitness ratio in original environment
			lin_norecalc_filename = run_directory + "/data/lineage_domEQU_norecalc.dat"
			lnfd = open(lin_norecalc_filename, "r")
			
			norecalc = lnfd.readlines()
			for line_num, line in enumerate(norecalc):
				line = line.strip()
				if len(line) == 0 or line[0] == "#":
					continue
				line = line.split()
				
				# We are looking specifically for the update when EQU fixed.
				# I am lazily looping through the lines.
				if int(line[7]) != fixed_EQU or int(line[7]) == -1:
					continue

				prev_line = norecalc[line_num - 1].strip().split()
				original_fixed_fitness = float(line[5]) / float(prev_line[5])
			
			lnfd.close()
		
	
	# add this data to tasks_info
	record.append(first_EQU)
	record.append(fixed_EQU)
	record.append(first_change)
	record.append(fixed_change)
	record.append(fixed_fitness)
	record.append(original_fixed_fitness)
	tasks_info.append(record)
	
	
# Output new file
output_filename = args[0] + "_taskinfo.dat"
ofd = open(output_filename, 'w')

ofd.write("# Avida - Aquisition of EQU task data\n")
ofd.write("# If an event never happens, the update is noted as -1; task change is 0.\n")
ofd.write("# 1: Seed\n")
ofd.write("# 2: Update at which EQU first appeared in population\n")

if options.use_lineage:
  ofd.write("# 3: Update at which EQU fixes in the final dominant's lineage\n")
else:
  ofd.write("# 3: Update at which EQU appears in population to stay until end-of-run\n")
  
ofd.write("# 4: Task change at first appearance of EQU (over population)\n")

if options.use_lineage:
  ofd.write("# 5: Task change at fixation of EQU (in the final dominant's lineage)\n")
else:
  ofd.write("# 5: Task change at fixation of EQU (over population)\n")
  
if options.use_lineage:
  ofd.write("# 6: Fitness ratio at fixation of EQU\n")
else:
 	ofd.write("# 6: N/A\n")
 	
if options.use_lineage:
	ofd.write("#7: Fitness ratio at fixation of EQU in original environment\n")
else:
	ofd.write("#7: N/A\n")
  
ofd.write("\n")

for record in tasks_info:
  ofd.write(str(record[0]) + " " + str(record[1]) + " " + str(record[2]) + " " + str(record[3]) + " " + str(record[4]) + " " + str(record[5]) + " " + str(record[6]) + "\n")
  

