# This script has two purposes:
# 	1. Collect the data regarding how many tasks are done per organism and per population over time.
#	   This is IO intensive.
#	2. Graph that data.  This is not so IO intensive.
#
# Therefore step 1 includes saving the collected information in a less IO intensive format.
# Actually doing step 1 is optional as long as this data exists; we may want to regenerate
# with step 2 many times without changing the backing data.  (i.e. different markers, colors, etc.)
# 
# Any number of experiments may be graphed at the same time.  Collected data for each run
# is stored in <experiment-name>_tasksper.dat, which has the format:
# <update> <median of tasks per organism> <median of tasks per population> <5th percentile per organism> <95 percentile per organism> <5th percentile per population> <95th percentile per population>
#

# Pretty pictures please
import matplotlib as mpl
mpl.use("TkAgg")

from argparse import ArgumentParser
import os.path
import numpy
import scipy.stats as ss
from pylab import *

# Constants
NUM_UPDATES = 100000
UPDATE_INTERVAL = 500


# Set up options and help
parser = ArgumentParser(description = "Collect and graph tasks per organism and per population over time.")

parser.add_argument("-r", "--repeat_collection", action = "store_true",
					default = False, help = "Read the original .spop files for all populations \
					over all updates, generating new data for tasks per organism and tasks per \
					population and recording them in <experiment_name>_tasksper.dat.  This is \
					applies to each experiment specified.  This is IO intensive and should be done \
					only when necessary; it happens automatically on a per-experiment basis when \
					a _tasksper.dat file for that experiment cannot be found. (default: %(default)s).")
parser.add_argument("-g", "--graph", action = "store_true", default = False,
					help = "Graph tasks per org and per population over time for each given \
					experiment (all on the same graph). (default: %(default)s)")
parser.add_argument("experiment_names", metavar = "E", nargs = "+",
					help = "An experiment name (e.g. unlimited-base) to be for which data will \
					be graphed or a _tasksper.dat file will be generated.")
					
# Of course, parsing takes place only if this file is called as a script; see the end of this file.

# Helper functions

# We only know the experiment names, not the seeds.
# We deduce which directories hold data of interest by
# examining all files/directories next to this script, and selecting:
# 1. Things where the experiment name matches one of ours.
# 	 a. We deduce this by stripping off the final underscore 
#       and what follows it.
# 2. Things that end in a seed -- so we know it's a run.
#    b. We simply look at what we stripped.  Is it underscore + number?
def find_run_directories(exp_name):	
	run_dirs = [dir for dir in os.listdir('./') 
	            if "_".join(dir.split("_")[:-1]) == exp_name 
	            and dir.split("_")[-1].isdigit()]
	            
	return run_dirs
	
	
# Collecting data and saving it in _tasksper.dat files.

# For each update timepoint,
#	for each run,
#		read the appropriate phenotypes file
#		gathering per-organism data and total per-population data
#	find the median, 5th, and 9th percentiles for per-o and per-p
#	write the update and the data to the _tasksper.dat file.
#
# input file format (<experiment_name>_phenotypes/phenotypes-<experiment_name>_seed-update.dat):
# <# orgs> <# genotypes> <avg. genome length> <avg. gest time> <viability> <binary task signature, space-separated>
#
# output file format (<experiment_name>_tasksper.dat):
# <update> <median of tasks per organism> <median of tasks per population> <5th percentile per organism> <95 percentile per organism> <5th percentile per population> <95th percentile per population>
def generate_tasksper_file(exp_name):

	# Write header on the _tasksper.dat file
	dat_fp = open(exp_name + "_tasksper.dat", "w")
	
	dat_fp.write("# 1. Update\n")
	dat_fp.write("# 2. Median number of distinct tasks performed per organism\n")
	dat_fp.write("# 3. Median number of distinct tasks performed per population\n")
	dat_fp.write("# 4. 5th percentile, per organism\n")
	dat_fp.write("# 5. 95th percentile, per organism\n")
	dat_fp.write("# 6. 5th percentile, per population\n")
	dat_fp.write("# 7. 95th percentile, per population\n")
	dat_fp.write("\n")
	
	dat_fp.close()
	
	# Phenotype files aren't in the run directories, but in <exp_name>_phenotypes/
	phen_dir = exp_name + "_phenotypes"
	
	# Derive seeds from experiment name
	run_dirs = find_run_directories(exp_name)
	seeds = [int(run_name.split("_")[-1]) for run_name in run_dirs]
	seeds.sort()
	
	for update in range(UPDATE_INTERVAL, NUM_UPDATES + 1, UPDATE_INTERVAL):
		org_tot_num_cpus = numpy.zeros(10, dtype=numpy.dtype(int))
		pop_tot_num_cpus  = numpy.zeros(10, dtype=numpy.dtype(int))

		for seed in seeds:
			# We get phenotype binary task signatures from the file,
			# but we must also collect a binary task signature for the entire population
			# as well as the population's total number of organisms
			pop_task_signature = [0] * 9
			pop_size = 0
			
			phen_fp = open(phen_dir + "/phenotypes-" + exp_name + "-" + str(seed) + "-" + str(update) + ".dat")
			
			for line in phen_fp:
				line = line.strip()
				if len(line) == 0 or line[0] == '#':  # skip empty lines and comments
					continue
				line = line.split()
				
				# How many tasks does this phenotype do?
				# Count the number of "1"s in the binary task signature
				phen_task_signature = [int(t) for t in line[5:]]
				num_tasks = sum(phen_task_signature)
				
				# Update the population's binary task signature
				pop_task_signature = [pop_task or phen_task for (pop_task, phen_task) in zip(pop_task_signature, phen_task_signature)]
				
				# How many organisms are in this phenotype?
				num_orgs = int(line[0])
				
				# Accumulate the number of organisms in the population
				pop_size += num_orgs

				# Data collected in bins, later exploded into a distribution to do stats
				org_tot_num_cpus[num_tasks] += num_orgs
				
			pop_tot_num_cpus[sum(pop_task_signature)] += pop_size
			
		# Now we've collected phenotype and population level data for this update
		# Explode it into a distribution, calculate the stats, 
		# and write a line to the _tasksper.dat file
		
		org_task_distribution = numpy.arange(10).repeat(org_tot_num_cpus)
		
		org_5th_percentile   = ss.scoreatpercentile(org_task_distribution, 5)
		org_median           = numpy.median(org_task_distribution)
		org_95th_percentile  = ss.scoreatpercentile(org_task_distribution, 95)
		
		pop_task_distribution = numpy.arange(10).repeat(pop_tot_num_cpus)
		
		pop_5th_percentile  = ss.scoreatpercentile(pop_task_distribution, 5)
		pop_median          = numpy.median(pop_task_distribution)
		pop_95th_percentile = ss.scoreatpercentile(pop_task_distribution, 95)
		
		# Write
		dat_fp = open(exp_name + "_tasksper.dat", "a")
		dat_fp.write(str(update) + " " + str(org_median) + " " + str(pop_median) \
		                         + " " + str(org_5th_percentile) + " " + str(org_95th_percentile) \
		                         + " " + str(pop_5th_percentile) +  " " + str(pop_95th_percentile) \
		                         + "\n")
		dat_fp.close()
	
	return
	
# Graphing!
def plot_tasksper(exp_names):

	for exp_name in exp_names:
		updates = []
		org_medians = []
		pop_medians = []
		
		# Read _tasksper.dat file, grabbing data
		fp = open(exp_name + "_tasksper.dat")
		
		for line in fp:
			line = line.strip()
			if len(line) == 0 or line[0] == '#':
				continue
			line = line.split()
			
			updates.append(int(line[0]))
			org_medians.append(float(line[1]))
			pop_medians.append(float(line[2]))
			
		fp.close()
		
		# Graph
		plot(updates, pop_medians, linewidth=2.0, alpha=0.7 ,label="per pop, "+exp_name)
		plot(updates, org_medians, linestyle='--', linewidth=3.0, alpha=0.7, label="per org, "+exp_name)
		
	legend(loc='lower center')
	title("Median # of Unique Tasks")
	show()
			
			

# When this file is run as a script:
if __name__ == "__main__":
	args = parser.parse_args()
	args = vars(args)	# transform Namespace into dictionary
	
	if args["repeat_collection"]:
		for experiment in args["experiment_names"]:
			generate_tasksper_file(experiment)
			
	if args["graph"]:
		for experiment in args["experiment_names"]:
			if not os.path.exists(experiment + "_tasksper.dat"):
				print "No _tasksper.dat file was found for " + experiment + "; generating."
				generate_tasksper_file(experiment)
		plot_tasksper(args["experiment_names"])