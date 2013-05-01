# Reads <experiment_name>_taskinfo.dat files, deduces correct population
# detail files to read, reads said files, and produces a file with lines
# in this format:
#
# seed phentropy_at_first_EQU phentropy_at_fixed_EQU phentropy_at_end
#
# phentropies are marked as -1 if EQU did not appear or did not fix
#
# Note that phenotypic entropy is not generated from the population at
# exactly the update of first_EQU and fixed_EQU, but the previous update
# divisible by 50.  If greater accuracy is desired, runs should be 
# resubmitted with an events file that dumps the population detail
# only at first_EQU and fixed_EQU.
#
# Bess L. Walker
# 8-20-09
# 9-15-09
# 7-11-11
#
# Note that I use the coinage "phentropy" because I am far too lazy to 
# write "phenotypic entropy".

import math
from argparse import ArgumentParser
import os.path

# Set up options
parser = ArgumentParser(description = "Given an experiment name, uses the associated _taskinfo.dat and phenotype files to generate a file describing the phenotypic entropy for each seed at the first appearance of EQU, the fixation of EQU, and the end of the run.")

parser.add_argument("experiment_name")
parser.add_argument("-f", "--filter-inviable", action = "store_true", 
                  dest = "filter_inviable", default = False, 
                  help = "include only viable genotypes in phenotypic entropy calculation (default: %(default)s)")
                  
args = parser.parse_args()
args = vars(args) # turn Namespace object into dictionary

experiment_name = args["experiment_name"]
filter_inviable = args["filter_inviable"]

phenotypes_dir = os.path.abspath(experiment_name + "_phenotypes")
taskinfo_filename = experiment_name + "_taskinfo.dat"

update_period = 50
num_updates = 100000

# Helper function that grabs the approximated phenotypes file, opens it
# read-only, and returns the file descriptor
# True update (truevalue) is approximated by the nearest previous update
# divisible by 50.  (If truevalue is divisible by 50, this is truevalue.)
#
# Caller is responsible for closing the file descriptor once it has
# been received.
def openapprox(seed, truevalue):
	
	difference = truevalue % 50
	approxvalue = truevalue - difference
	
	approx_filename = phenotypes_dir + "/phenotypes-" + experiment_name + "-" + str(seed) + "-" + str(approxvalue) + ".dat"
	
	return open(approx_filename, "r")
	

# Helper function which calculates phentropy given a phenotypefile
# as file descriptor; does not close file
def phentropy(pfd):
	org_counts = []
	
	for line in pfd:
		line = line.strip()
		
		if len(line) == 0 or line[0] == "#":
			continue
			
		line = line.split()
			
		# if using only viable phenotypes, filter for them
		if not filter_inviable or int(line[4]) != 0:
			org_counts.append(int(line[0]))
	
	total_orgs = float(sum(org_counts))
	org_fractions = [count / total_orgs for count in org_counts]
	entropy = -1 * sum([frac * math.log10(frac) for frac in org_fractions])
	
	return entropy
	

# Everything else is here:

phentropy_info = []

# Read tasks file and phenotypes files; calculate phentropies

tfd = open(taskinfo_filename, "r")

for line in tfd:
	line = line.strip()
	
	if len(line) == 0 or line[0] == "#":
		continue
		
	line = line.split()
	
	seed = int(line[0])
	first_EQU = int(line[1])
	fixed_EQU = int(line[2])
	
	if first_EQU != -1:
		first_fd = openapprox(seed, first_EQU)
		first_phentropy = phentropy(first_fd)
		first_fd.close()
	else:
		first_phentropy = -1
		
	if fixed_EQU != -1:
		fixed_fd = openapprox(seed, fixed_EQU)
		fixed_phentropy = phentropy(fixed_fd)
		fixed_fd.close()
	else:
		fixed_phentropy = -1
		
	final_fd = openapprox(seed, num_updates)	
	final_phentropy = phentropy(final_fd)
	final_fd.close()
	
	phentropy_info.append([seed, first_phentropy, fixed_phentropy, final_phentropy])
	
tfd.close()

# Write

# use a different filename for inviable-filtered phentropy
if filter_inviable:
	out_filename = experiment_name + "_phentropyinfo_viable.dat"
else:
	out_filename = experiment_name + "_phentropyinfo.dat"
	
ofd = open(out_filename, "w")

ofd.write("# Avida - Aquisition of EQU phenotypic entropy data\n")
ofd.write("# If the event never happened, phenotypic entropy is noted as -1.\n")
ofd.write("# 1. Seed\n")
ofd.write("# 2. Phenotypic entropy at which EQU first appeared in population\n")
ofd.write("# 3. Phenotypic entropy at which EQU appears to stay until end-of-run\n")
ofd.write("# 4. Phenotpyic entropy at end of run\n")
ofd.write("\n")

for record in phentropy_info:
	ofd.write(str(record[0]) + " ")
	ofd.write(str(record[1]) + " ")
	ofd.write(str(record[2]) + " ")
	ofd.write(str(record[3]) + "\n")

ofd.close()
		