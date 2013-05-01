# Calculates and displays statistics relating to the ecovo experiments:
#
# As appropriate, minimum, median, and maxium are reported; for statistical
# tests p-values and any other relevant figures are reported.
#
# NUMBER OF RUNS:
# - ever aquiring EQU
# - fixing EQU
# - having EQU at end-of-run (should be the same as number fixing)
# 
# TIME:
# - to first aquire EQU
# - to fix EQU
#
# TASK TRADEOFFS:
# - # in the population when it first acquires EQU
# - # in the dominant EQU lineage when it fixes EQU
#
# PHENOTYPIC ENTROPY:
# - at first aquistion of EQU
# - at fixation of EQU
# - at end of run
# (Phenotypic entropy is reported with and without including inviable 
# phenotypes.)
#
# When we speak of "fixation" of EQU we mean the first update in which the
# population performed EQU and did not lose EQU in any subsequent update; this
# is opposed to the first aquisition of EQU, which may be subsequently lost.
# 
#
#
#
# WHAT TO RUN:
#   To see p-values based on filtering out runs which did not achieve, run:
#   > do_reports(numerical_filter)
#
#   To see p-values based on all runs (non-achieve = inf), run:
#   > do_reports(lambda L: L)
#
# The variables exp_1_names and exp_2_names should be defined 
# (i.e. "unlimited-base" and "limited-base") if the script is not called 
# from the command line.
#
#
#
# Written with Python 2.7.1; requires 2.6 or greater
# Requires SciPy 0.9.0 or greater for Fisher's Exact Test
# BLW 8-31-09
# BLW 9-15-09
# BLW 9-22-09
# BLW 7-12-11, 7-13-11, 7-14-11

import scipy
import numpy
import scipy.stats.stats as sss
import scipy.stats as ss
from argparse import ArgumentParser
import os
import os.path
import gzip

# SET UP OPTIONS AND HELP
usage = "usage: %prog [options]"
parser = ArgumentParser(description = "Calculate and report various numbers and statistics for the experiments compared.  By default, runs which did not achieve a particular milestone (e.g. first EQU, fixed EQU) are filtered out of the calculations for that milestone.")

parser.add_argument("-u", "--unfiltered", action = "store_true", 
                  default = False, help = "Include all runs in all statistics, using a value of infinity for any milestone which was not achieved (default: %(default)s)")
parser.add_argument("-e1", "--experiment_1_names", metavar = "NAME", required = True, nargs = "*", help = "The name(s) of the first experiment, i.e. unlimited-base.  One name is typical.")
parser.add_argument("-e2", "--experiment_2_names", metavar = "NAME", required = True, nargs = "*", help = "The name(s) of the second experiment, i.e. limited-base.  One name is typical.")
                  
# Parsing takes place only if this is called as a script; see end of file                  


# DEFINE CONSTANTS

FINAL_UPDATE = 100000
RUNS_PER_EXPERIMENT = 200

LABEL = "LABEL"
FIRST = "FIRST"
FIXED = "FIXED"
FINAL = "FINAL"
UPDATES = "UPDATES"
TRADEOFF = "TRADEOFF"
FITNESS = "FITNESS RATIO"
ORIG_FITNESS = "ORIGINAL FITNESS RATIO"
ENTROPY = "ENTROPY"
V_ENTROPY = "V_ENTROPY"
TASKS_PER_ORG = "TASKS PER ORGANISM"
LOSS_NO_GAIN = "TASK LOSS WITHOUT GAIN"
LNG_BEFORE_EQU = "TASK LOSS WITHOUT GAIN DIRECTLY PRECEEDING EQU GAIN"
EQU_PER_1000 = "NEW EQU EVOLVED PER 1000 UPDATES"
EQU_PER_1000_CONSERVATIVE = "NEW EQU EVOLVED PER 1000 UPDATES AFTER FIRST EQU"
TOT_GENOTYPES_EQU = "TOTAL NUMBER OF GENOTYPES (POP. EVOLVED EQU)"
TOT_GENOTYPES_NOEQU = "TOTAL NUMBER OF GENOTYPES (POP. DID NOT EVOLVE EQU)" 
DATA = "DATA"
COUNT = "COUNT"
MIN = "MIN"
PER_5 = "FIFTH PERCENTILE"
MEDIAN = "median "
PER_95 = "NINETY-FIFTH PERCENTILE"
MAX = "MAX"

HIGH_UPDATE = scipy.inf
HIGH_TRADEOFF = scipy.inf
HIGH_PHENTROPY = scipy.inf
HIGH_FITNESS = scipy.inf

HIGH_LIST = [HIGH_UPDATE, HIGH_TRADEOFF, HIGH_PHENTROPY, HIGH_FITNESS]

# DEFINE HELPER FUNCTIONS

# Creates an empty dictionary to hold all the data we are going to collect
def make_ecoevo_stats_dict():
  return { FIRST : { UPDATES : {}, ENTROPY : {}, V_ENTROPY: {}}, 
           FIXED : { UPDATES : {}, TRADEOFF : {}, FITNESS : {}, 
                     ORIG_FITNESS: {}, ENTROPY : {}, V_ENTROPY: {}}, 
           FINAL : { ENTROPY : {}, V_ENTROPY: {}, TASKS_PER_ORG: {}, 
                     LOSS_NO_GAIN: {}, LNG_BEFORE_EQU: {}, 
                     EQU_PER_1000: {}, EQU_PER_1000_CONSERVATIVE: {},
                     TOT_GENOTYPES_EQU: {}, TOT_GENOTYPES_NOEQU: {}}}
                    
# Finds the maximum actual number in a list that also contains HIGH values
def numerical_max(L):
  return numpy.max([x for x in L if (x not in HIGH_LIST)])
  
# Finds the fifth percentile of the actual numbers in a list that also 
# contains HIGH values
def numerical_5_percentile(L):
    return ss.scoreatpercentile([x for x in L if (x not in HIGH_LIST)], 5)

# Finds the median actual number in a list that also contains HIGH values
def numerical_median(L):
  return numpy.median([x for x in L if (x not in HIGH_LIST)])
  
# Finds the 95th percentile of the actual numbers in a list that also
# contains HIGH values
def numerical_95_percentile(L):
    return ss.scoreatpercentile([x for x in L if (x not in HIGH_LIST)], 95)
  
# Returns a list of the non-HIGH numbers in a list
def numerical_filter(L):
  return [x for x in L if (x != HIGH_UPDATE and x != HIGH_PHENTROPY)]
  

# READ FILES
# Update information and phenotypic entropy information are read in
# and translated to formats best suited for further analysis (i.e.
# -1 meaning "never happened" becomes inf).

# This function reads an <experiment_name>_taskinfo.dat file that 
# contains update numbers for first_EQU and fixed_EQU, as well as 
# the number of tasks gained/lost at each milestone.
# 
# Arguments: experiment_names (list of experiment names that belong to
# the same experiment)
#
# Returns: first_EQU_updates, fixed_EQU_updates, first_EQU_tradeoffs, 
#          fixed_EQU_tradeoffs (lists of updates and tradeoffs
#          from this experiment)
def read_taskfile(experiment_names):
  first_EQU_updates = []
  fixed_EQU_updates = []
  fixed_EQU_tradeoffs = []
  fixed_EQU_fitness = []
  fixed_EQU_original_fitness = []
  
  for exp in experiment_names:
    efd = open(exp + "_taskinfo.dat", "r")
        
    for line in efd:
      line = line.strip()
      if len(line) == 0 or line[0] == "#":
        continue
      line = line.split()
      
      # read and store first EQU update #
      if int(line[1]) < 0:
        first_EQU_updates.append(HIGH_UPDATE)
      else:
        first_EQU_updates.append(int(line[1]))
        
      # read and store fixed EQU update #
      if int(line[2]) < 0:
        fixed_EQU_updates.append(HIGH_UPDATE)
      else:
        fixed_EQU_updates.append(int(line[2]))
        
      # don't read tradeoff at first EQU; it's meaningless
      # (Well, it has meaning at a population level, but 
      # it's not very interesting.)
      
      # read and store fixed EQU tradeoff
      if int(line[2]) < 0:
        fixed_EQU_tradeoffs.append(HIGH_TRADEOFF)
      else:
        fixed_EQU_tradeoffs.append(int(line[4]))
        
      # read and store fixed EQU fitness ratio
      if float(line[5]) < 0:
        fixed_EQU_fitness.append(HIGH_FITNESS)
      else:
        fixed_EQU_fitness.append(float(line[5]))
        
      # read and store fixed EQU fitness ratio from original environment
      if float(line[6]) < 0:
        fixed_EQU_original_fitness.append(HIGH_FITNESS)
      else:
        fixed_EQU_original_fitness.append(float(line[6]))
    
    efd.close()
  
  return first_EQU_updates, fixed_EQU_updates, fixed_EQU_tradeoffs, fixed_EQU_fitness, fixed_EQU_original_fitness
  
# This function reads an <experiment_name>_phentropyinfo.dat or 
# <experiment_name>_phentropyinfo_viable.dat file that contains 
# phenotypic entropy at first_EQU, fixed_EQU, and end of run.
#
# Arguments: 
#   - experiment_names (list of experiment names that belong to the
# same experiment)
#   - viable_only (True = read _viable.dat, False = read _phentropyinfo.dat)
#
# Returns: first_EQU_phentropies, fixed_EQU_phentropies, final_EQU_phentropies
# (lists of phenotypic entropies from this experiment)
def read_phentropyfile(experiment_names, viable_only):
  first_EQU_phentropies = []
  fixed_EQU_phentropies = []
  final_EQU_phentropies = []
  
  for exp in experiment_names:
    if viable_only:
      efd = open(exp + "_phentropyinfo_viable.dat", "r")
    else:
      efd = open(exp + "_phentropyinfo.dat", "r")
      
    for line in efd:
      line = line.strip()
      if len(line) == 0 or line[0] == "#":
        continue
      line = line.split()
      
      # read and store phentropy at first EQU
      if float(line[1]) < 0:
        first_EQU_phentropies.append(HIGH_PHENTROPY)
      else:
        first_EQU_phentropies.append(float(line[1]))
        
      # read and store phentropy at fixed EQU
      if float(line[2]) < 0:
        fixed_EQU_phentropies.append(HIGH_PHENTROPY)
      else:
        fixed_EQU_phentropies.append(float(line[2]))
        
      # read and store phentropy at end of run;
      # this always exists
      final_EQU_phentropies.append(float(line[3]))
        
    efd.close()
  
  return first_EQU_phentropies, fixed_EQU_phentropies, final_EQU_phentropies

# Reads phenotype files to determine the distribution of # tasks performed.  
def read_phenotypefiles(experiment_names):
    tot_num_cpus = numpy.zeros(10, dtype=numpy.dtype(int))
    
    # Phenotypes are in <experiment_name>_phenotypes/
    # They are named phenotypes-<experiment_name>-seed-update.dat
    # We want all the update 100000 files.
    for exp in experiment_names:
        phenfiles = [file for file in os.listdir(exp + "_phenotypes") if file.split("-")[-1] == "100000.dat"]
        
        for filename in phenfiles:
            pfd = open(exp + "_phenotypes/" + filename, "r")
            
            for line in pfd:
                line = line.strip()
                if len(line) == 0 or line[0] == '#':
                    continue
                line = line.split()
                
                # read number of cpus for this phenotype
                num_cpus = int(line[0])
                
                # read task signature for this phenotype
                tasks = [int(t) for t in line[5:]]
                
                # count number of tasks
                num_tasks = sum(tasks)
                
                # store this data
                tot_num_cpus[num_tasks] += num_cpus
            
            pfd.close()
            
    # Now we have all the data, explode it into a distribution:
    num_task_distribution = numpy.arange(10).repeat(tot_num_cpus)
    
    return num_task_distribution
    
# Reads lineage files and counts the number of task-loss-without-task-gain
# mutations, and the number of these that occur directly before EQU acquisition.
# No entry appears for empty lineages (where the population didn't achieve EQU).
def read_lineagefiles(experiment_names):
    loss_no_gains = []
    lng_before_EQU = []
    lng_beneficials = []
    lng_bene_before_EQU = []
    
    # We only know the experiment names, not the seeds.
    # We deduce which directories hold data of interest by
    # examining all files/directories next to this script, and selecting:
    # 1. Things where the experiment name matches one of ours.
    #    a. We deduce this by stripping off the final underscore 
    #       and what follows it.
    # 2. Things that end in a seed -- so we know it's a run.
    #    b. We simply look at what we stripped.  Is it underscore + number?
    
    run_dirs = [dir for dir in os.listdir('./') 
                if "_".join(dir.split("_")[:-1]) in experiment_names 
                and dir.split("_")[-1].isdigit()]
    
    for dir in run_dirs:
        if not os.path.exists(dir + "/data/lineage_domEQU.dat"):
            continue
        
        # Grab the data from the non-recalculated lineage
        # These manipulations mean that indexing into norecalc_lines
        # can be done using the line_count that goes through the legitimate
        # lines of the recalculated lineage (below).
        norecalc_lfd = open(dir + "/data/lineage_domEQU_norecalc.dat", "r")
        norecalc_lines = norecalc_lfd.readlines()
        norecalc_lines = [line.strip() for line in norecalc_lines]
        norecalc_lines = [line for line in norecalc_lines if len(line) != 0 and line[0] != "#"]
        norecalc_lines = [line.split() for line in norecalc_lines]
            
        # Actually work with the data from the recalculated lineage 
        lfd = open(dir + "/data/lineage_domEQU.dat", "r")
        line_count = -1
        grandparent_tasks = []
        num_lng = 0
        num_lng_EQU = 0
        num_bene_lng = 0
        num_bene_lng_EQU = 0
        
        for line in lfd:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            else:
                line_count += 1
            
            line = line.split()
            
            # Is this an empty lineage?  Check that depth matches the line number.
            if int(line[0]) != line_count:
                break

            # The ancestor does no tasks:
            if line_count == 0:
                grandparent_tasks   = [0] * 9
                grandparent_fitness = float(norecalc_lines[line_count][5])
                grandparent_gest    = float(norecalc_lines[line_count][6])
            # We can't possibly lose tasks with the first mutation, since we had 0 anyway.
            elif line_count == 1:
                parent_tasks    = [int(t) for t in line[8:]]
                parent_fitness  = float(norecalc_lines[line_count][5])
                parent_gest     = float(norecalc_lines[line_count][6])
            else:
                # We compare parent to grandparent for loss without gain
                # If we see loss without gain, we look at the current task to see if it is EQU
                
                current_tasks = [int(t) for t in line[8:]]
                
                # Gains show up as 1, losses as -1, no change as 0
                task_change = [p - g for p, g in zip(parent_tasks, grandparent_tasks)]
                
                # No gain, and also loss?
                if 1 not in task_change and -1 in task_change:
                
                    # Do we think this is probably a non-real loss,
                    # i.e. caused by stochastic performance?
                    #
                    # That is, we consider a mutation to likely be a real loss if:
                    #   a) it was detrimental
                    #   b) it was neutral or beneficial, and gestation time decreased
                    #
                    # This is not perfect, but we are limited by the data.
                    if parent_fitness < grandparent_fitness or parent_gest < grandparent_gest:
                        num_lng += 1
                        if parent_fitness >= grandparent_fitness:
                            num_bene_lng += 1
                        # Gain EQU directly after this mutation?
                        if current_tasks[-1] == 1 and parent_tasks[-1] != 1:
                            num_lng_EQU += 1
                            if parent_fitness >= grandparent_fitness:
                                num_bene_lng_EQU += 1
                                
                grandparent_tasks = parent_tasks[:]
                parent_tasks = current_tasks[:]
                grandparent_fitness = parent_fitness
                parent_fitness = float(norecalc_lines[line_count][5])
                grandparent_gest = parent_gest
                parent_gest = float(norecalc_lines[line_count][6])
                        
        loss_no_gains.append(num_lng)
        lng_before_EQU.append(num_lng_EQU)
        lng_beneficials.append(num_bene_lng)
        lng_bene_before_EQU.append(num_bene_lng_EQU)
                        
        lfd.close()
        
    print "There are " + str(len([x for x in lng_beneficials if x != 0])) + " lineages with probably-real beneficial LnGs."
    print "Of these, " + str(len([x for x in lng_bene_before_EQU if x != 0])) + " had a beneficial LnG directly before EQU."
    print

    return loss_no_gains, lng_before_EQU

# Reads newtasks.dat files and calculates the number of new EQU evolutions
# per 1000 updates.  Further calculates the number of new EQU evolutions per
# 1000 updates only counting updates after the first EQU.
# Concurrently, reads totals.dat and records the total number of genotypes
# in the run, building separate lists for runs that did get EQU and runs that didn't.
def read_newtaskfiles(experiment_names):
    per_1000 = []
    per_1000_conservative = []
    tot_genotypes_EQU = []
    tot_genotypes_noEQU = []
    
    # First we deduce the run directories as in read_lineagefiles():
    run_dirs = [dir for dir in os.listdir('./')
                if "_".join(dir.split("_")[:-1]) in experiment_names
                and dir.split("_")[-1].isdigit()]
    
    # Then we read the files, if they are there.
    # If they are not there, we warn.
    for dir in run_dirs:
        if not os.path.exists(dir + "/data/newtasks.dat.gz"):
#           print "Warning: Unable to find newtasks.dat.gz file for " + dir
            continue
        
        nfd = gzip.open(dir + "/data/newtasks.dat.gz", "r")
        first_EQU_update = -1
        EQU_sum = 0
        
        for line in nfd:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            line = line.split()
            
            # Any new EQU here?
            if int(line[9]) > 0:
                EQU_sum += int(line[9])
                # First EQU encountered?
                if first_EQU_update == -1:
                    first_EQU_update = int(line[0])
                    
        # Calculate figures for this run, add to lists
        if first_EQU_update == -1:
            per_1000.append(0)  
            per_1000_conservative.append(0)  # these don't count, but they are filtered out later in populate_EQUper1000_stats()
        else:
            per_1000.append(float(EQU_sum) / float(FINAL_UPDATE) * 1000)
            updates_after_EQU = FINAL_UPDATE - first_EQU_update
            per_1000_conservative.append(float(EQU_sum - 1) / float(updates_after_EQU) * 1000)
            
        nfd.close()
            
        # Now read totals.dat to find the total number of genotypes in this run.
        if not os.path.exists(dir + "/data/totals.dat.gz"):
            print "Warning: Unable to find totals.dat.gz file for " + dir
            continue
            
        tfd = gzip.open(dir + "/data/totals.dat.gz", "r")
        
        # We only care about the last line, at update 100000
        for line in tfd:
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                continue
            line = line.split()
            if line[0] != "100000":
                continue
                
            # Okay, how many genotypes had been produced by update 100000?
            tot_genotypes = int(line[4])
            
            # If the population never got EQU, the first update with EQU is -1.
            if first_EQU_update == -1:
                tot_genotypes_noEQU.append(tot_genotypes)
            else:
                tot_genotypes_EQU.append(tot_genotypes)
        
        tfd.close()

    return per_1000, per_1000_conservative, tot_genotypes_EQU, tot_genotypes_noEQU
    
    
  
# Reads all necessary files and stores lists to DATA elements
def populate(dict, experiment_names):
    # Even if there's more than one name for this experiment, we have to choose.
    # So we arbitrarily choose the first one, then prettify it.
    dict[LABEL] = experiment_names[0].replace("-base", "").capitalize()
    
    
    dict[FIRST][UPDATES][DATA],  \
    dict[FIXED][UPDATES][DATA],  \
    dict[FIXED][TRADEOFF][DATA], \
    dict[FIXED][FITNESS][DATA],  \
    dict[FIXED][ORIG_FITNESS][DATA] = read_taskfile(experiment_names)
    
    dict[FIRST][ENTROPY][DATA], \
    dict[FIXED][ENTROPY][DATA], \
    dict[FINAL][ENTROPY][DATA] = read_phentropyfile(experiment_names, False)
    
    dict[FIRST][V_ENTROPY][DATA], \
    dict[FIXED][V_ENTROPY][DATA], \
    dict[FINAL][V_ENTROPY][DATA] = read_phentropyfile(experiment_names, True)
    
    dict[FINAL][TASKS_PER_ORG][DATA] = read_phenotypefiles(experiment_names)

    dict[FINAL][LOSS_NO_GAIN][DATA], \
    dict[FINAL][LNG_BEFORE_EQU][DATA] = read_lineagefiles(experiment_names)
    
    dict[FINAL][EQU_PER_1000][DATA], \
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][DATA], \
    dict[FINAL][TOT_GENOTYPES_EQU][DATA], \
    dict[FINAL][TOT_GENOTYPES_NOEQU][DATA] = read_newtaskfiles(experiment_names)

    
    
# NUMBER OF RUNS
# Numbers of runs for each EQU-point are counted and reported;
# they are compared with scipy's two-tailed Fisher's Exact Test

# Counts non-HIGH_UPDATE values in a list of update values
def count_achieved(update_list):
  return len([x for x in update_list if x != HIGH_UPDATE])
  
# Does all the counts
def populate_counts(dict):
  dict[FIRST][UPDATES][COUNT] = count_achieved(dict[FIRST][UPDATES][DATA])
  dict[FIXED][UPDATES][COUNT] = count_achieved(dict[FIXED][UPDATES][DATA])
  
# Compare with Fisher's Exact Test (two-tailed)
# No filtering necessary for this one.
def do_count_comparisons(exp_1, exp_2):
    first_table = [
                  [exp_1[FIRST][UPDATES][COUNT], 
                   RUNS_PER_EXPERIMENT - exp_1[FIRST][UPDATES][COUNT]],
                  [exp_2[FIRST][UPDATES][COUNT],
                   RUNS_PER_EXPERIMENT - exp_2[FIRST][UPDATES][COUNT]]
                  ]
    fixed_table = [
                  [exp_1[FIXED][UPDATES][COUNT], 
                   RUNS_PER_EXPERIMENT - exp_1[FIXED][UPDATES][COUNT]],
                  [exp_2[FIXED][UPDATES][COUNT],
                   RUNS_PER_EXPERIMENT - exp_2[FIXED][UPDATES][COUNT]]
                  ] 
                  
    first = ss.fisher_exact(first_table)
    
    if fixed_table[0] == fixed_table[1]: # There's no difference
        fixed = (1, 1)
    else:
        fixed = ss.fisher_exact(fixed_table)
    
    return first, fixed

# Report
def report_counts(exp_1, exp_2, label_width, first_fet, fixed_fet):
    print
    print "NUMBER OF RUNS (achieving this milestone):"
    print
    print "{:<{width}} | First EQU | Fixed EQU".format('', width = label_width)
    print "{:{fill}<{width}}-|-----------|----------" \
          .format('', fill = "-", width = label_width)
          
    # Experiment 1 numbers      
    print "{:<{width}} | {:^9} | {:^9}" \
          .format(exp_1[LABEL], exp_1[FIRST][UPDATES][COUNT], exp_1[FIXED][UPDATES][COUNT], width = label_width)
    
    # Experiment 2 numbers
    print "{:<{width}} | {:^9} | {:^9}" \
          .format(exp_2[LABEL], exp_2[FIRST][UPDATES][COUNT], exp_2[FIXED][UPDATES][COUNT], width = label_width)
          
    print "{:{fill}<{width}}-|-----------|----------" \
          .format('', fill = "-", width = label_width)
    
    # p-value
    print "{:<{width}} | {:9.3e} | {:9.3e}" \
          .format(" p-value", first_fet[1], fixed_fet[1], width = label_width)
          
    print
    print "(Uses Fisher's exact test, two-tailed)"
    print
    


# TIME
# Minimum, median, and maximum # updates to each EQU-point are 
# reported; the distributions are compared via Mann-Whitney test
# and p-values are reported.

def populate_time_stats(dict):
    first_list = dict[FIRST][UPDATES][DATA]
    fixed_list = dict[FIXED][UPDATES][DATA]
    
    # If there are no first EQUs in any population...
    if len(numerical_filter(first_list)) == 0:
        dict[FIRST][UPDATES][MIN]    = 0
        dict[FIRST][UPDATES][PER_5]  = 0
        dict[FIRST][UPDATES][MEDIAN] = 0
        dict[FIRST][UPDATES][PER_95] = 0
        dict[FIRST][UPDATES][MAX]    = 0
    else:
        dict[FIRST][UPDATES][MIN]    = min(first_list)
        dict[FIRST][UPDATES][PER_5] = numerical_5_percentile(first_list)
        dict[FIRST][UPDATES][MEDIAN] = numerical_median(first_list)
        dict[FIRST][UPDATES][PER_95] = numerical_95_percentile(first_list)
        dict[FIRST][UPDATES][MAX]    = numerical_max(first_list)
        
    # If there are no fixed EQUs in any population...
    if len(numerical_filter(fixed_list)) == 0:
        dict[FIXED][UPDATES][MIN]    = 0
        dict[FIXED][UPDATES][PER_5]  = 0
        dict[FIXED][UPDATES][MEDIAN] = 0
        dict[FIXED][UPDATES][PER_95] = 0
        dict[FIXED][UPDATES][MAX]    = 0
    else:
        dict[FIXED][UPDATES][MIN]    = min(fixed_list)
        dict[FIXED][UPDATES][PER_5] = numerical_5_percentile(fixed_list)
        dict[FIXED][UPDATES][MEDIAN] = numerical_median(fixed_list)
        dict[FIXED][UPDATES][PER_95] = numerical_95_percentile(fixed_list)
        dict[FIXED][UPDATES][MAX]    = numerical_max(fixed_list)

# Compare with Mann-Whitney U
def do_time_comparisons(filter_func, exp_1, exp_2):
  first = sss.mannwhitneyu(filter_func(exp_1[FIRST][UPDATES][DATA]), 
                           filter_func(exp_2[FIRST][UPDATES][DATA]))
  
  # if no population fixed EQU
  if len(numerical_filter(exp_1[FIXED][UPDATES][DATA])) == 0 and \
     len(numerical_filter(exp_2[FIXED][UPDATES][DATA])) == 0: 
    fixed = (1, 1)
  else:
    fixed = sss.mannwhitneyu(filter_func(exp_1[FIXED][UPDATES][DATA]), 
                           filter_func(exp_2[FIXED][UPDATES][DATA]))
  return first, fixed                                    
                                      
                                    
# Report times
def report_times(exp_1, exp_2, label_width, first_mwu, fixed_mwu):
  print
  print "TIME (in updates to achieve this milestone):"
  print
  print "{:<{width}} | First EQU | Fixed EQU".format('', width = label_width)
  print "{:{fill}<{width}}-|-----------|----------" \
        .format('', fill = '-', width = label_width)
  
  # Experiment 1
  print "{:<{width}} |           |          " \
        .format(exp_1[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("min    ", exp_1[FIRST][UPDATES][MIN],     \
                exp_1[FIXED][UPDATES][MIN],     width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("5 %ile ", exp_1[FIRST][UPDATES][PER_5], \
                exp_1[FIXED][UPDATES][PER_5],  width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("median ", exp_1[FIRST][UPDATES][MEDIAN], \
                 exp_1[FIXED][UPDATES][MEDIAN], width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("95 %ile", exp_1[FIRST][UPDATES][PER_95], \
                exp_1[FIXED][UPDATES][PER_95],  width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("max    ", exp_1[FIRST][UPDATES][MAX],    \
                exp_1[FIXED][UPDATES][MAX],     width = label_width)
  
  # Experiment 2
  print "{:<{width}} |           |          " \
        .format(exp_2[LABEL], width = label_width)
  
  print "{:>{width}} | {:9d} | {:9d}"    \
        .format("min    ", exp_2[FIRST][UPDATES][MIN],     \
                exp_2[FIXED][UPDATES][MIN],     width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("5 %ile ", exp_2[FIRST][UPDATES][PER_5], \
                exp_2[FIXED][UPDATES][PER_5],  width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}" \
        .format("median ", exp_2[FIRST][UPDATES][MEDIAN], \
                 exp_2[FIXED][UPDATES][MEDIAN], width = label_width)
  print "{:>{width}} | {:9.1f} | {:9.1f}"     \
        .format("95 %ile", exp_2[FIRST][UPDATES][PER_95], \
                exp_2[FIXED][UPDATES][PER_95],  width = label_width)
  print "{:>{width}} | {:9d} | {:9d}"     \
        .format("max    ", exp_2[FIRST][UPDATES][MAX],    \
                exp_2[FIXED][UPDATES][MAX],     width = label_width)
  
                                       
  print "{:{fill}<{width}}-|-----------|----------" \
        .format('', fill = '-', width = label_width)
  
  # p-value
  print "{:>{width}} | {:9.3e} | {:9.3e}" \
        .format(" p-value", first_mwu[1], fixed_mwu[1], \
                width = label_width)
  print
  print "(Uses Mann-Whitney U test, one-tailed)"
  print



# TASKS PER ORGANISM
# The number of tasks per organism at the end of the run.
#
# We compare the distributions with the Mann-Whitney U-test and
# report p-values, median, max, and min.

def populate_tasksper_stats(dict):
    task_list = dict[FINAL][TASKS_PER_ORG][DATA]
    
    # We do NOT use the numerical_x() functions defined earlier for
    # two reasons:
    # 1. The number of tasks ranges from 0 to 9 inclusive.  There are
    #    by definition no HIGH values.
    # 2. Running a list comprehension on something the size of
    #    task_list (roughly 640,000 entries) takes a rather
    #    noticeable amount of time.
    #
    # We DO make very sure to use numpy-based functions.
    
    dict[FINAL][TASKS_PER_ORG][MIN]    = numpy.min(task_list)
    dict[FINAL][TASKS_PER_ORG][PER_5] = ss.scoreatpercentile(task_list, 5)
    dict[FINAL][TASKS_PER_ORG][MEDIAN] = numpy.median(task_list)
    dict[FINAL][TASKS_PER_ORG][PER_95] = ss.scoreatpercentile(task_list, 95)
    dict[FINAL][TASKS_PER_ORG][MAX]    = numpy.max(task_list)
        
# Compare with Mann-Whitney U
def do_tasksper_comparisons(exp_1, exp_2):
    final = sss.mannwhitneyu(exp_1[FINAL][TASKS_PER_ORG][DATA], 
                             exp_2[FINAL][TASKS_PER_ORG][DATA])
    
    return final

# Report tasks per organism
def report_tasksperorg(exp_1, exp_2, label_width, final_mwu):
    print
    print "TASKS PER ORGANISM:"
    print
    print "{:<{width}} | End of Run".format('', width = label_width)
    print "{:{fill}<{width}}-|----------" \
                .format('', fill = "-", width = label_width)
    
    # Experiment 1
    print "{:<{width}} |          ".format(exp_1[LABEL], width = label_width)
    
    print "{:>{width}} | {:9d}"   \
          .format("min    ", exp_1[FINAL][TASKS_PER_ORG][MIN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("5 %ile ", exp_1[FINAL][TASKS_PER_ORG][PER_5], 
                  width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("median ", exp_1[FINAL][TASKS_PER_ORG][MEDIAN], 
                  width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("95 %ile", exp_1[FINAL][TASKS_PER_ORG][PER_95], 
                  width = label_width)
    print "{:>{width}} | {:9d}" \
          .format("max    ", exp_1[FINAL][TASKS_PER_ORG][MAX], width = label_width)
    
    # Experiment 2
    print "{:<{width}} |          ".format(exp_2[LABEL], width = label_width)
    
    print "{:>{width}} | {:9d}" \
          .format("min    ", exp_2[FINAL][TASKS_PER_ORG][MIN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("5 %ile ", exp_2[FINAL][TASKS_PER_ORG][PER_5], 
                  width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("median ", exp_2[FINAL][TASKS_PER_ORG][MEDIAN], 
                  width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("95 %ile", exp_2[FINAL][TASKS_PER_ORG][PER_95], 
                  width = label_width)
    print "{:>{width}} | {:9d}" \
          .format("max    ", exp_2[FINAL][TASKS_PER_ORG][MAX], width = label_width)
    
    print "{:{fill}<{width}}-|----------" \
                .format('', fill = "-", width = label_width)
    
    # p-value
    print "{:<{width}} | {:e}".format(" p-value", final_mwu[1], width = label_width)
    
    print
    print "(Uses Mann-Whitney U test, one-tailed)"
    print

    
# TASK TRADEOFFS
# Number of tasks gained (+) or lost (-) at each milestone.
# The data for FIXED is the only meaningful data here, as it comes from 
# a single lineage (the final dominant EQU) rather than the whole
# population -- a population which may be full of all sorts of 
# genotypes gaining and losing tasks willy-nilly.
#
# We compare the distribution with the Mann-Whitney U-test and report
# p-values.

def populate_tradeoff_stats(dict):
  fixed_list = dict[FIXED][TRADEOFF][DATA]
  
  # Are all tradeoffs infinite?  i.e. there were no runs with EQU
  if len(numerical_filter(fixed_list)) == 0:
        dict[FIXED][TRADEOFF][MIN]    = 0
        dict[FIXED][TRADEOFF][PER_5] = 0
        dict[FIXED][TRADEOFF][MEDIAN] = 0
        dict[FIXED][TRADEOFF][PER_95] = 0
        dict[FIXED][TRADEOFF][MAX]    = 0
        return
    

  dict[FIXED][TRADEOFF][MIN]    = min(fixed_list)
  dict[FIXED][TRADEOFF][PER_5] = numerical_5_percentile(fixed_list)
  dict[FIXED][TRADEOFF][MEDIAN] = numerical_median(fixed_list)
  dict[FIXED][TRADEOFF][PER_95] = numerical_95_percentile(fixed_list)
  dict[FIXED][TRADEOFF][MAX]    = numerical_max(fixed_list)
  
# Compare with Mann-Whitney U
def do_tradeoff_comparisons(filter_func, exp_1, exp_2):
    
    # If no population fixed EQU...
    if len(filter_func(exp_1[FIXED][TRADEOFF][DATA])) == 0 and \
       len(filter_func(exp_2[FIXED][TRADEOFF][DATA])) == 0:
        fixed = (1, 1)
    else:
        fixed = sss.mannwhitneyu(filter_func(exp_1[FIXED][TRADEOFF][DATA]),
                                                     filter_func(exp_2[FIXED][TRADEOFF][DATA]))
                                                     
    return fixed
  
# Report tradeoffs
def report_tradeoffs(exp_1, exp_2, label_width, fixed_mwu):
    print
    print "TRADEOFFS (number tasks gained/lost at this milestone):"
    print
    print "{:<{width}} | Fixed EQU".format('', width = label_width)
    print "{:{fill}<{width}}-|----------" \
                .format('', fill = "-", width = label_width)
    
    # Experiment 1
    print "{:<{width}} |          ".format(exp_1[LABEL], width = label_width)
    
    print "{:>{width}} | {:9d}"   \
          .format("min    ", exp_1[FIXED][TRADEOFF][MIN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("5 %ile ", exp_1[FIXED][TRADEOFF][PER_5], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("median ", exp_1[FIXED][TRADEOFF][MEDIAN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("95 %ile", exp_1[FIXED][TRADEOFF][PER_95], width = label_width)
    print "{:>{width}} | {:9d}"   \
          .format("max    ", exp_1[FIXED][TRADEOFF][MAX], width = label_width)
    
    # Experiment 2
    print "{:<{width}} |          ".format(exp_2[LABEL], width = label_width)
    
    print "{:>{width}} | {:9d}"   \
          .format("min    ", exp_2[FIXED][TRADEOFF][MIN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("5 %ile ", exp_2[FIXED][TRADEOFF][PER_5], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("median ", exp_2[FIXED][TRADEOFF][MEDIAN], width = label_width)
    print "{:>{width}} | {:9.1f}" \
          .format("95 %ile", exp_2[FIXED][TRADEOFF][PER_95], width = label_width)
    print "{:>{width}} | {:9d}"   \
          .format("max    ", exp_2[FIXED][TRADEOFF][MAX], width = label_width)
    
    print "{:{fill}<{width}}-|----------" \
                .format('', fill = "-", width = label_width)
    
    # p-value
    print "{:<{width}} | {:9.3e}".format(" p-value", fixed_mwu[1], width = label_width)
    
    print
    print "(Uses Mann-Whitney U test, one-tailed)"
    print


# FITNESS RATIO
# Fitness ratio of parent to child when child is the genotype which 
# acquired the fixing EQU, where fitness is evaluated always in an
# unlimited resources environment (the test CPU).  So this is a 
# "what-would-have-been" for the limited-resource data.
#
# Again, since this is based on lineage, we have data only for 
# fixation.
#
# We compare the distributions using a Mann-Whitney U test.

def populate_fitness_stats(dict):
  fixed_list = dict[FIXED][FITNESS][DATA]
  fixed_orig_list = dict[FIXED][ORIG_FITNESS][DATA]
  
  if len(numerical_filter(fixed_list)) == 0:
    return
  
  dict[FIXED][FITNESS][MIN]    = min(fixed_list)
  dict[FIXED][FITNESS][PER_5]  = numerical_5_percentile(fixed_list)
  dict[FIXED][FITNESS][MEDIAN] = numerical_median(fixed_list)
  dict[FIXED][FITNESS][PER_95] = numerical_95_percentile(fixed_list)
  dict[FIXED][FITNESS][MAX]    = numerical_max(fixed_list)
  
  dict[FIXED][ORIG_FITNESS][MIN]    = min(fixed_orig_list)
  dict[FIXED][ORIG_FITNESS][PER_5]  = numerical_5_percentile(fixed_orig_list)
  dict[FIXED][ORIG_FITNESS][MEDIAN] = numerical_median(fixed_orig_list)
  dict[FIXED][ORIG_FITNESS][PER_95] = numerical_95_percentile(fixed_orig_list)
  dict[FIXED][ORIG_FITNESS][MAX]    = numerical_max(fixed_orig_list)
  
# Compare with Mann-Whitney U
def do_fitness_comparisons(filter_func, exp_1, exp_2):
  if len(filter_func(exp_1[FIXED][FITNESS][DATA])) == 0 or len(filter_func(exp_2[FIXED][FITNESS][DATA])) == 0:
    return (-1, -1)
    
  fixed = sss.mannwhitneyu(filter_func(exp_1[FIXED][FITNESS][DATA]),
                           filter_func(exp_2[FIXED][FITNESS][DATA]))
                           
  fixed_orig = sss.mannwhitneyu(filter_func(exp_1[FIXED][ORIG_FITNESS][DATA]),
                                filter_func(exp_2[FIXED][ORIG_FITNESS][DATA]))
                           
  return fixed, fixed_orig
  
# Report fitness ratios
def report_fitness_ratios(exp_1, exp_2, label_width, fixed_mwu):
  if fixed_mwu == -1:
    print "Fitness ratios not available; data does not derive from a lineage"
    return

  print
  print "FITNESS RATIO (evaluated in unlimited envrionment):"
  print
  print "{:>{width}} | Fixed EQU".format('', width = label_width)
  print "{:{fill}>{width}}-|----------" \
        .format('', fill = '-', width = label_width)
        
  # Experiment 1
  print "{:<{width}} |          ".format(exp_1[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.1f}" \
        .format("min    ", exp_1[FIXED][FITNESS][MIN],    width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("5 %ile ", exp_1[FIXED][FITNESS][PER_5], width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("median ", exp_1[FIXED][FITNESS][MEDIAN], width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("95 %ile", exp_1[FIXED][FITNESS][PER_95], width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("max    ", exp_1[FIXED][FITNESS][MAX],    width = label_width)
  
  # Experiment 2
  print "{:<{width}} |          ".format(exp_2[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.1f}"   \
        .format("min    ", exp_2[FIXED][FITNESS][MIN],    width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("5 %ile ", exp_2[FIXED][FITNESS][PER_5], width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("median ", exp_2[FIXED][FITNESS][MEDIAN], width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("95 %ile", exp_2[FIXED][FITNESS][PER_95], width = label_width)
  print "{:>{width}} | {:9.1f}"   \
        .format("max    ", exp_2[FIXED][FITNESS][MAX],    width = label_width)
  
  print "{:{fill}>{width}}-|----------" \
        .format('', fill = '-', width = label_width)
  
  # p-value -- multiplied by 2 to get two-tailed value
  print "{:<{width}} | {:9.3e}" \
        .format(" p-value", fixed_mwu[1] * 2.0, width = label_width)
  
  print
  print "(Uses Mann-Whitney U test, two-tailed)"
  print
  
# Report fitness ratios in original environment
def report_orig_fitness_ratios(exp_1, exp_2, label_width, fixed_orig_mwu):
  if fixed_orig_mwu == -1:
    print "Fitness ratios not available; data does not derive from a lineage"
    return

  print
  print "FITNESS RATIO (evaluated in original environment):"
  print
  print "{:>{width}} | Fixed EQU".format('', width = label_width)
  print "{:{fill}>{width}}-|----------" \
        .format('', fill = '-', width = label_width)
        
  # Experiment 1
  print "{:<{width}} |          ".format(exp_1[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.1f}" \
        .format("min    ", exp_1[FIXED][ORIG_FITNESS][MIN],
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("5 %ile ", exp_1[FIXED][ORIG_FITNESS][PER_5], 
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("median ", exp_1[FIXED][ORIG_FITNESS][MEDIAN], 
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("95 %ile", exp_1[FIXED][ORIG_FITNESS][PER_95], 
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("max    ", exp_1[FIXED][ORIG_FITNESS][MAX],    
                width = label_width)
  
  # Experiment 2
  print "{:<{width}} |          ".format(exp_2[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.1f}"   \
        .format("min    ", exp_2[FIXED][ORIG_FITNESS][MIN],    
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("5 %ile ", exp_2[FIXED][ORIG_FITNESS][PER_5], 
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("median ", exp_2[FIXED][ORIG_FITNESS][MEDIAN], 
                width = label_width)
  print "{:>{width}} | {:9.1f}" \
        .format("95 %ile", exp_2[FIXED][ORIG_FITNESS][PER_95], 
                width = label_width)
  print "{:>{width}} | {:9.1f}"   \
        .format("max    ", exp_2[FIXED][ORIG_FITNESS][MAX],    
                 width = label_width)
  
  print "{:{fill}>{width}}-|----------" \
        .format('', fill = '-', width = label_width)
  
  # p-value -- multiplied by 2.0 to get two-tailed value
  print "{:<{width}} | {:9.3e}" \
        .format(" p-value", fixed_orig_mwu[1] * 2.0, width = label_width)
  
  print
  print "(Uses Mann-Whitney U test, two-tailed)"
  print
  
  


# TASK LOSS WITHOUT TASK GAIN (on dominant EQU-achieving lineage)
# We count the number of lineages on which this specific type of detrimental
# mutation appears.
#
# We compare with Fisher's Exact Test, and report the p-value.

def populate_lossnogain_stats(dict):
    loss_list = dict[FINAL][LOSS_NO_GAIN][DATA]
    equ_list  = dict[FINAL][LNG_BEFORE_EQU][DATA]
    
    dict[FINAL][LOSS_NO_GAIN][COUNT] = len([x for x in loss_list if x != 0])
    dict[FINAL][LNG_BEFORE_EQU][COUNT] = len([x for x in equ_list if x != 0])
        
def do_lossnogain_comparisons(exp_1, exp_2):
    exp_1_runs_with_EQU = exp_1[FIXED][UPDATES][COUNT]
    exp_2_runs_with_EQU = exp_2[FIXED][UPDATES][COUNT]

    if exp_1_runs_with_EQU == 0 and exp_2_runs_with_EQU == 0:
        final_lng = (1, 1)
        final_equ = (1, 1)
        return final_lng, final_equ

    final_lng_table = [
                       [exp_1[FINAL][LOSS_NO_GAIN][COUNT],
                        exp_1_runs_with_EQU - exp_1[FINAL][LOSS_NO_GAIN][COUNT]],
                       [exp_2[FINAL][LOSS_NO_GAIN][COUNT],
                        exp_2_runs_with_EQU - exp_2[FINAL][LOSS_NO_GAIN][COUNT]]
                      ]
                      
    final_equ_table = [
                       [exp_1[FINAL][LNG_BEFORE_EQU][COUNT],
                        exp_1[FINAL][LOSS_NO_GAIN][COUNT] - exp_1[FINAL][LNG_BEFORE_EQU][COUNT]],
                       [exp_2[FINAL][LNG_BEFORE_EQU][COUNT],
                        exp_2[FINAL][LOSS_NO_GAIN][COUNT] - exp_2[FINAL][LNG_BEFORE_EQU][COUNT]]
                      ]
        
    if final_lng_table[0] == final_lng_table[1]:
        final_lng = (1, 1)
    else:
        final_lng = ss.fisher_exact(final_lng_table)      

    if final_equ_table[0] == final_equ_table[1]:
        final_equ = (1, 1)
    else:
        final_equ = ss.fisher_exact(final_equ_table)
    
    return final_lng, final_equ 
    
def report_lossnogain(exp_1, exp_2, label_width, final_lng_fet, final_equ_fet):
    print
    print "TASK LOSS WITHOUT TASK GAIN:"
    print
    print "{:<{width}} | Number of Runs".format("", width = label_width)
    print "{:{fill}<{width}}-|---------------" \
          .format('', fill = "-", width = label_width)
          
    # Experiment 1 numbers      
    print "{:<{width}} | {:^14} " \
          .format(exp_1[LABEL], exp_1[FINAL][LOSS_NO_GAIN][COUNT], 
                  width = label_width)
    
    # Experiment 2 numbers
    print "{:<{width}} | {:^14}" \
          .format(exp_2[LABEL], exp_2[FINAL][LOSS_NO_GAIN][COUNT],
                  width = label_width)
          
    print "{:{fill}<{width}}-|----------------" \
          .format('', fill = "-", width = label_width)
    
    # p-value
    print "{:<{width}} | {:14.3e}" \
          .format(" p-value", final_lng_fet[1], width = label_width)
          
    print
    print "(Uses Fisher's exact test, two-tailed)"
    print
    
    print
    print "EQU MUTATIONS DIRECTLY FOLLOWING TASK LOSS WITHOUT TASK GAIN:"
    print
    print "{:<{width}} | Number of Runs".format("", width = label_width)
    print "{:{fill}<{width}}-|---------------" \
          .format('', fill = "-", width = label_width)
          
    # Experiment 1 numbers      
    print "{:<{width}} | {:^14} " \
          .format(exp_1[LABEL], exp_1[FINAL][LNG_BEFORE_EQU][COUNT], 
                  width = label_width)
    
    # Experiment 2 numbers
    print "{:<{width}} | {:^14}" \
          .format(exp_2[LABEL], exp_2[FINAL][LNG_BEFORE_EQU][COUNT],
                  width = label_width)
          
    print "{:{fill}<{width}}-|----------------" \
          .format('', fill = "-", width = label_width)
    
    # p-value
    print "{:<{width}} | {:14.3e}" \
          .format(" p-value", final_equ_fet[1], width = label_width)
          
    print
    print "(Uses Fisher's exact test, two-tailed)"
    print

    
    
# NUMBER OF EQU EVOLUTIONS PER 1000 UPDATES
# We look at how often EQU newly evolves.
# The distributions are compared with the Mann-Whitney U test.

def populate_EQUper1000_stats(dict):
    final_list = dict[FINAL][EQU_PER_1000][DATA]    
    final_c_list = dict[FINAL][EQU_PER_1000_CONSERVATIVE][DATA]
    
    if len(final_list) == 0 or len(final_c_list) == 0:
        print "Statistics for EQU evolutions per 1000 are not available."
        return
    
    # Conservative measure does not include non-EQU-evolving runs
    final_c_list = [x for x in final_c_list if x != 0]
    
    dict[FINAL][EQU_PER_1000][MIN]    = min(final_list)
    dict[FINAL][EQU_PER_1000][PER_5]  = numerical_5_percentile(final_list)
    dict[FINAL][EQU_PER_1000][MEDIAN] = numerical_median(final_list)
    dict[FINAL][EQU_PER_1000][PER_95] = numerical_95_percentile(final_list)
    dict[FINAL][EQU_PER_1000][MAX]    = numerical_max(final_list)
    
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][MIN]    = min(final_c_list)
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][PER_5]  = numerical_5_percentile(final_c_list)
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][MEDIAN] = numerical_median(final_c_list)
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][PER_95] = numerical_95_percentile(final_c_list)
    dict[FINAL][EQU_PER_1000_CONSERVATIVE][MAX]    = numerical_max(final_c_list)
    
def do_EQUper1000_comparisons(filter_func, exp_1, exp_2):
    # If the newtasks.dat was not available:
    if len(exp_1[FINAL][EQU_PER_1000][DATA]) == 0 or \
       len(exp_2[FINAL][EQU_PER_1000][DATA]) == 0 or \
       len(exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][DATA]) == 0 or \
       len(exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][DATA]) == 0:
        return (-1, -1), (-1, -1)

    final = sss.mannwhitneyu(filter_func(exp_1[FINAL][EQU_PER_1000][DATA]),
                             filter_func(exp_2[FINAL][EQU_PER_1000][DATA]))
                             
    final_conservative = sss.mannwhitneyu(filter_func(exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][DATA]),
                     filter_func(exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][DATA]))
                     
    return final, final_conservative
    
def report_EQUper1000(exp1, exp2, label_width, final_mwu, final_c_mwu):
  if final_mwu == -1:
    print "Frequency of EQU evolution not available; new task data not available for this data set."
    return

  print
  print "EQU EVOLUTIONS PER 1000 UPDATES:"
  print
  print "{:>{width}} | {:^9} | {:12}".format('', "EQU/1000", "CONSERVATIVE", width = label_width)
  print "{:{fill}>{width}}-|-----------|--------------" \
        .format('', fill = '-', width = label_width)
  
  # If these stats were not calculated
  if final_mwu == (-1, -1) or final_c_mwu == (-1, -1):
    print "Data for these statistics was not calculated, and therefore they cannot be displayed."
    return
        
  # Experiment 1
  print "{:<{width}} |          ".format(exp_1[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("min    ", exp_1[FINAL][EQU_PER_1000][MIN],
                exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][MIN],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("5 %ile ", exp_1[FINAL][EQU_PER_1000][PER_5], 
                exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][PER_5],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("median ", exp_1[FINAL][EQU_PER_1000][MEDIAN],
                exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][MEDIAN],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("95 %ile", exp_1[FINAL][EQU_PER_1000][PER_95],
                exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][PER_95],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("max    ", exp_1[FINAL][EQU_PER_1000][MAX],
                exp_1[FINAL][EQU_PER_1000_CONSERVATIVE][MAX],
                width = label_width)
  
  # Experiment 2
  print "{:<{width}} |          ".format(exp_2[LABEL], width = label_width)
  
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("min    ", exp_2[FINAL][EQU_PER_1000][MIN],
                exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][MIN],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("5 %ile ", exp_2[FINAL][EQU_PER_1000][PER_5], 
                exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][PER_5],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("median ", exp_2[FINAL][EQU_PER_1000][MEDIAN],
                exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][MEDIAN],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("95 %ile", exp_2[FINAL][EQU_PER_1000][PER_95],
                exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][PER_95],
                width = label_width)
  print "{:>{width}} | {:9.2f} | {:12.3f}" \
        .format("max    ", exp_2[FINAL][EQU_PER_1000][MAX],
                exp_2[FINAL][EQU_PER_1000_CONSERVATIVE][MAX],
                width = label_width)
  
  print "{:{fill}>{width}}-|----------|--------------" \
        .format('', fill = '-', width = label_width)
  
  # p-values -- multiplied by 2 to get the two-tailed p-values
  print "{:<{width}} | {:9.3e} | {:14.3e}" \
        .format(" p-value", final_mwu[1] * 2.0, final_c_mwu[1] * 2.0, width = label_width)
  
  print
  print "(Uses Mann-Whitney U test, two-tailed)"
  print
  
  

# NUMBER OF TOTAL GENOTYPES IN EQU- AND NON-EQU-EVOLVING RUNS
# We look at how many total genotypes were generated.
# The distributions are compared WITHIN treatments with the Mann-Whitney U test.
# (This data is only generated for the EQU-resampled treatments.)

def populate_totgenotypes_stats(dict):
    final_EQU_list = dict[FINAL][TOT_GENOTYPES_EQU][DATA]
    final_noEQU_list = dict[FINAL][TOT_GENOTYPES_NOEQU][DATA]
    
    if len(final_EQU_list) == 0 and len(final_noEQU_list) == 0:
        print "Statistics for total number of genotypes are not available."
        return
        
    for (label, genotype_list) in [(TOT_GENOTYPES_EQU, final_EQU_list), \
                                   (TOT_GENOTYPES_NOEQU, final_noEQU_list)]:                 
        if len(genotype_list) != 0:                           
            dict[FINAL][label][MIN]    = min(genotype_list)
            dict[FINAL][label][PER_5]  = numerical_5_percentile(genotype_list)
            dict[FINAL][label][MEDIAN] = numerical_median(genotype_list)
            dict[FINAL][label][PER_95] = numerical_95_percentile(genotype_list)
            dict[FINAL][label][MAX]    = numerical_max(genotype_list)
        else:
            dict[FINAL][label][MIN]    = numpy.nan
            dict[FINAL][label][PER_5]  = numpy.nan
            dict[FINAL][label][MEDIAN] = numpy.nan
            dict[FINAL][label][PER_95] = numpy.nan
            dict[FINAL][label][MAX]    = numpy.nan           

def do_totgenotypes_comparisons(filter_func, exp_1, exp_2):
    final_counts = []
    for exp in [exp_1, exp_2]:
        # If the statistics are not available:
        if len(exp[FINAL][TOT_GENOTYPES_EQU][DATA]) == 0 and \
           len(exp[FINAL][TOT_GENOTYPES_NOEQU][DATA]) == 0:
            final_counts.append((-1, -1))
        else:
            final = sss.mannwhitneyu(filter_func(exp[FINAL][TOT_GENOTYPES_EQU][DATA]),
                                     filter_func(exp[FINAL][TOT_GENOTYPES_NOEQU][DATA]))
            final_counts.append(final)
            
    return final_counts[0], final_counts[1]
    
def report_totgenotypes(exp_1, exp_2, label_width, final_exp1_mwu, final_exp2_mwu):
    print
    print "TOTAL GENOTYPES OVER THE RUN:"
    print
    print "{:>{width}} | {:^11} | {:^12}".format('', "EQU", "NO EQU", width = label_width)
    print "{:{fill}>{width}}-|-------------|--------------" \
    .format('', fill = '-', width = label_width)

    for (exp, final_mwu) in [(exp_1, final_exp1_mwu), (exp_2, final_exp2_mwu)]:            
        # If these stats were not calculated:
        if final_mwu == (-1, -1):
            print "Total number of genotypes per run was not calculated for", exp[LABEL] + "."
            continue

        print "{:<{width}} |          ".format(exp[LABEL], width = label_width)

        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("min    ", exp[FINAL][TOT_GENOTYPES_EQU][MIN],
                      exp[FINAL][TOT_GENOTYPES_NOEQU][MIN],
                      width = label_width)
        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("5 %ile ", exp[FINAL][TOT_GENOTYPES_EQU][PER_5], 
                      exp[FINAL][TOT_GENOTYPES_NOEQU][PER_5],
                      width = label_width)
        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("median ", exp[FINAL][TOT_GENOTYPES_EQU][MEDIAN],
                      exp[FINAL][TOT_GENOTYPES_NOEQU][MEDIAN],
                      width = label_width)
        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("95 %ile", exp[FINAL][TOT_GENOTYPES_EQU][PER_95],
                      exp[FINAL][TOT_GENOTYPES_NOEQU][PER_95],
                      width = label_width)
        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("max    ", exp[FINAL][TOT_GENOTYPES_EQU][MAX],
                      exp[FINAL][TOT_GENOTYPES_NOEQU][MAX],
                      width = label_width)
        # p-value WITHIN treatment, multiplied by 2 to get the two-tailed value
        print
        print "{:>{width}} | {:11.2e} |" \
              .format ("p-value", final_mwu[1] * 2.0, width = label_width)
        print "{:>{width}} | {:11.2f} | {:12.3f}" \
              .format("min    ", len(exp[FINAL][TOT_GENOTYPES_EQU][DATA]),
                      len(exp[FINAL][TOT_GENOTYPES_NOEQU][DATA]),
                      width = label_width)
        print


    print
    print "(Uses Mann-Whitney U test, two-tailed)"
    print

        

# PHENOTYPIC ENTROPY:
# Minimum, median, and maximum entropy at each EQU-point are reported;
# the distributions are compared via Mann-Whitney test and p-values are
# reported.

def populate_entropy_stats(dict):
  first_list = dict[FIRST][ENTROPY][DATA]
  first_v_list = dict[FIRST][V_ENTROPY][DATA]
  fixed_list = dict[FIXED][ENTROPY][DATA]
  fixed_v_list = dict[FIXED][V_ENTROPY][DATA]

  dict[FIRST][ENTROPY][MIN]    = min(first_list)
  dict[FIRST][ENTROPY][PER_5]  = numerical_5_percentile(first_list)
  dict[FIRST][ENTROPY][MEDIAN] = numerical_median(first_list)
  dict[FIRST][ENTROPY][PER_95] = numerical_95_percentile(first_list)
  dict[FIRST][ENTROPY][MAX]    = numerical_max(first_list)
  
  dict[FIRST][V_ENTROPY][MIN]    = min(first_v_list)
  dict[FIRST][V_ENTROPY][PER_5]  = numerical_5_percentile(first_v_list)
  dict[FIRST][V_ENTROPY][MEDIAN] = numerical_median(first_v_list)
  dict[FIRST][V_ENTROPY][PER_95] = numerical_95_percentile(first_v_list)
  dict[FIRST][V_ENTROPY][MAX]    = numerical_max(first_v_list)
  
  # If no populations fixed EQU...
  if len(numerical_filter(fixed_list)) == 0:
        dict[FIXED][ENTROPY][MIN]    = 0
        dict[FIXED][ENTROPY][PER_5]  = 0
        dict[FIXED][ENTROPY][MEDIAN] = 0
        dict[FIXED][ENTROPY][PER_95] = 0
        dict[FIXED][ENTROPY][MAX]    = 0
        
        dict[FIXED][V_ENTROPY][MIN]    = 0
        dict[FIXED][V_ENTROPY][PER_5]  = 0
        dict[FIXED][V_ENTROPY][MEDIAN] = 0
        dict[FIXED][V_ENTROPY][PER_95] = 0
        dict[FIXED][V_ENTROPY][MAX]    = 0
  else: 
        dict[FIXED][ENTROPY][MIN]    = min(fixed_list)
        dict[FIXED][ENTROPY][PER_5]  = numerical_5_percentile(fixed_list)
        dict[FIXED][ENTROPY][MEDIAN] = numerical_median(fixed_list)
        dict[FIXED][ENTROPY][PER_95] = numerical_95_percentile(fixed_list)
        dict[FIXED][ENTROPY][MAX]    = numerical_max(fixed_list)
        
        dict[FIXED][V_ENTROPY][MIN]    = min(fixed_v_list)
        dict[FIXED][V_ENTROPY][PER_5]  = numerical_5_percentile(fixed_v_list)
        dict[FIXED][V_ENTROPY][MEDIAN] = numerical_median(fixed_v_list)
        dict[FIXED][V_ENTROPY][PER_95] = numerical_95_percentile(fixed_v_list)
        dict[FIXED][V_ENTROPY][MAX]    = numerical_max(fixed_v_list)
  
  populate_end_entropy_stats(dict)
  
def populate_end_entropy_stats(dict):
  final_list = dict[FINAL][ENTROPY][DATA]
  final_v_list = dict[FINAL][V_ENTROPY][DATA]
  
  dict[FINAL][ENTROPY][MIN]    = min(final_list)
  dict[FINAL][ENTROPY][PER_5]  = numerical_5_percentile(final_list)
  dict[FINAL][ENTROPY][MEDIAN] = numerical_median(final_list)
  dict[FINAL][ENTROPY][PER_95] = numerical_95_percentile(final_list)
  dict[FINAL][ENTROPY][MAX]    = numerical_max(final_list)
  
  dict[FINAL][V_ENTROPY][MIN]    = min(final_v_list)
  dict[FINAL][V_ENTROPY][PER_5]  = numerical_5_percentile(final_v_list)
  dict[FINAL][V_ENTROPY][MEDIAN] = numerical_median(final_v_list)
  dict[FINAL][V_ENTROPY][PER_95] = numerical_95_percentile(final_v_list)
  dict[FINAL][V_ENTROPY][MAX]    = numerical_max(final_v_list)
  
# Compare entropies using Mann-Whitney U test
def do_entropy_comparisons(entropy_type, filter_func, exp_1, exp_2):
  first = sss.mannwhitneyu(filter_func(exp_1[FIRST][entropy_type][DATA]), 
                           filter_func(exp_2[FIRST][entropy_type][DATA]))
  
  # if no populations fixed EQU:
  if len(filter_func(exp_1[FIXED][entropy_type][DATA])) == 0 and \
     len(filter_func(exp_2[FIXED][entropy_type][DATA])) == 0:
    fixed = (1, 1)
  else:
    fixed = sss.mannwhitneyu(filter_func(exp_1[FIXED][entropy_type][DATA]), 
                           filter_func(exp_2[FIXED][entropy_type][DATA]))

  final = sss.mannwhitneyu(filter_func(exp_1[FINAL][entropy_type][DATA]),
                           filter_func(exp_2[FINAL][entropy_type][DATA]))
                                          
  return first, fixed, final

# Report
def report_phentropy(label, key, exp_1, exp_2, label_width, first_mwu, fixed_mwu, final_mwu):
  print
  print "PHENOTYPIC ENTROPY (when achieving this milestone):"
  print
  print "{:<{width}} | First EQU | Fixed EQU | End Of Run" \
        .format(' ' + label, width = label_width)
  print "{:{fill}<{width}}-|-----------|-----------|-----------" \
        .format('', fill = '-', width = label_width)
  
  # Experiment 1
  print "{:<{width}} |           |           |           " \
        .format(exp_1[LABEL], width = label_width)
        
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("min    ", exp_1[FIRST][key][MIN], exp_1[FIXED][key][MIN], 
                exp_1[FINAL][key][MIN], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("5 %ile ", exp_1[FIRST][key][PER_5], exp_1[FIXED][key][PER_5],
                exp_1[FINAL][key][PER_5], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("median ", exp_1[FIRST][key][MEDIAN], exp_1[FIXED][key][MEDIAN], 
                exp_1[FINAL][key][MEDIAN], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("95 %ile", exp_1[FIRST][key][PER_95], 
                exp_1[FIXED][key][PER_95], exp_1[FINAL][key][PER_95], 
                width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("max    ", exp_1[FIRST][key][MAX], exp_1[FIXED][key][MAX], 
                exp_1[FINAL][key][MAX], width = label_width)
  
  # Experiment 2
  print "{:<{width}} |           |           |           " \
        .format(exp_2[LABEL], width = label_width)
        
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("min    ", exp_2[FIRST][key][MIN], exp_2[FIXED][key][MIN], 
                exp_2[FINAL][key][MIN], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("5 %ile ", exp_2[FIRST][key][PER_5], exp_2[FIXED][key][PER_5],
                exp_2[FINAL][key][PER_5], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("median ", exp_2[FIRST][key][MEDIAN], exp_2[FIXED][key][MEDIAN], 
                exp_2[FINAL][key][MEDIAN], width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("95 %ile", exp_2[FIRST][key][PER_95], 
                exp_2[FIXED][key][PER_95], exp_2[FINAL][key][PER_95], 
                width = label_width)
  print "{:>{width}} | {:.7f} | {:.7f} | {:.7f}" \
        .format("max    ", exp_2[FIRST][key][MAX], exp_2[FIXED][key][MAX], 
                exp_2[FINAL][key][MAX], width = label_width)
  
  print "{:{fill}<{width}}-|-----------|-----------|-----------" \
        .format('', fill = '-', width = label_width)
  
  # p-value -- multiplied by 2 to get two-tailed p-values
  print "{:<{width}} | {:.3e} | {:.3e} | {:.3e}".format(" p-value", first_mwu[1], fixed_mwu[1] * 2.0, final_mwu[1] * 2.0, width = label_width)  
  
  print
  print "(Uses Mann-Whitney U test, two-tailed)"
  print
              

# ALL STATS
# Functions that bring everything together

# This function does all data population of the dictionaries
def do_setup(exp_1_names, exp_2_names):
    # Set up dictionaries
    exp_1 = make_ecoevo_stats_dict()
    exp_2 = make_ecoevo_stats_dict()
    
    # Populate with data from files
    populate(exp_1, exp_1_names)
    populate(exp_2, exp_2_names)
    
    # Populate with EQU-achieving counts
    populate_counts(exp_1)
    populate_counts(exp_2)
    
    # Populate with task tradeoff stats
    populate_tradeoff_stats(exp_1)
    populate_tradeoff_stats(exp_2)
    
    # Populate with tasks per organism stats
    populate_tasksper_stats(exp_1)
    populate_tasksper_stats(exp_2)
    
    # Populate with fitness ratio at EQU fixation stats
    populate_fitness_stats(exp_1)
    populate_fitness_stats(exp_2)
    
    # Populate with stats about task gain without task loss
    populate_lossnogain_stats(exp_1)
    populate_lossnogain_stats(exp_2)
    
    # Populate with time-to-EQU stats
    populate_time_stats(exp_1)
    populate_time_stats(exp_2)
    
    # Populate with how often EQU evolves
    populate_EQUper1000_stats(exp_1)
    populate_EQUper1000_stats(exp_2)
    
    # Populate with how many total genotypes per run there are
    populate_totgenotypes_stats(exp_1)
    populate_totgenotypes_stats(exp_2)
    
    # Populate with phenotypic entropy stats
    populate_entropy_stats(exp_1)
    populate_entropy_stats(exp_2)
    
    # Compare whether phenotypic entropy with inviable genotypes and entropy
    # without inviable genotypes is actually different
#   exp_1_first_viable_mwu = sss.mannwhitneyu(exp_1[FIRST][ENTROPY][DATA], exp_1[FIRST][V_ENTROPY][DATA])
#   exp_1_fixed_viable_mwu = sss.mannwhitneyu(exp_1[FIXED][ENTROPY][DATA], exp_1[FIXED][V_ENTROPY][DATA])
#   exp_1_final_viable_mwu = sss.mannwhitneyu(exp_1[FINAL][ENTROPY][DATA], exp_1[FINAL][V_ENTROPY][DATA])
#   
#   exp_2_first_viable_mwu = sss.mannwhitneyu(exp_2[FIRST][ENTROPY][DATA],
#   exp_2[FIRST][V_ENTROPY][DATA])
#   exp_2_fixed_viable_mwu = sss.mannwhitneyu(exp_2[FIXED][ENTROPY][DATA],
#   exp_2[FIXED][V_ENTROPY][DATA])
#   exp_2_final_viable_mwu = sss.mannwhitneyu(exp_2[FINAL][ENTROPY][DATA],
#   exp_2[FINAL][V_ENTROPY][DATA])
    
    label_width = max(len(exp_1[LABEL]), len(exp_2[LABEL]))
    
    return exp_1, exp_2, label_width

def do_reports(filter_func, exp_1, exp_2, label_width):
    first_count_fet, fixed_count_fet = do_count_comparisons(exp_1, exp_2)
    report_counts(exp_1, exp_2, label_width, first_count_fet, fixed_count_fet)
    
    first_update_mwu, fixed_update_mwu = do_time_comparisons(filter_func, exp_1, exp_2)
    report_times(exp_1, exp_2, label_width, first_update_mwu, fixed_update_mwu)
    
    final_tasksper_mwu = do_tasksper_comparisons(exp_1, exp_2)
    report_tasksperorg(exp_1, exp_2, label_width, final_tasksper_mwu)
    
    fixed_tradeoff_mwu = do_tradeoff_comparisons(filter_func, exp_1, exp_2)
    report_tradeoffs(exp_1, exp_2, label_width, fixed_update_mwu)
    
    fixed_fitness_mwu, fixed_orig_fitness_mwu = do_fitness_comparisons(filter_func, exp_1, exp_2)
    report_fitness_ratios(exp_1, exp_2, label_width, fixed_fitness_mwu)
    report_orig_fitness_ratios(exp_1, exp_2, label_width, fixed_orig_fitness_mwu)
    
    final_lossnogain_fet, final_lossnogainequ_fet = do_lossnogain_comparisons(exp_1, exp_2)
    report_lossnogain(exp_1, exp_2, label_width, final_lossnogain_fet, final_lossnogainequ_fet)
    
    final_EQUper1000_mwu, final_EQUper1000conserv_mwu = do_EQUper1000_comparisons(filter_func, exp_1, exp_2)
    report_EQUper1000(exp_1, exp_2, label_width, final_EQUper1000_mwu, final_EQUper1000conserv_mwu)
    
    final_totgenotypesexp1_mwu, final_totgenotypesexp2_mwu = do_totgenotypes_comparisons(filter_func, exp_1, exp_2)
    report_totgenotypes(exp_1, exp_2, label_width, final_totgenotypesexp1_mwu, final_totgenotypesexp2_mwu)
    
    first_entropy_mwu, fixed_entropy_mwu, final_entropy_mwu = do_entropy_comparisons(ENTROPY, filter_func, exp_1, exp_2)
    report_phentropy("all orgs", ENTROPY, exp_1, exp_2, label_width, first_entropy_mwu, fixed_entropy_mwu, final_entropy_mwu)
                                     
    first_v_entropy_mwu, fixed_v_entropy_mwu, final_v_entropy_mwu = do_entropy_comparisons(V_ENTROPY, filter_func, exp_1, exp_2)
    report_phentropy("viable  ", V_ENTROPY, exp_1, exp_2, label_width, first_v_entropy_mwu, fixed_v_entropy_mwu, final_v_entropy_mwu)

# And in the end...
if __name__ == "__main__":
    args = parser.parse_args()
    args = vars(args) # transform Namespace into dictionary
    exp_1_names = args["experiment_1_names"]
    exp_2_names = args["experiment_2_names"]
    
    exp_1, exp_2, label_width = do_setup(exp_1_names, exp_2_names)
    
    if not args["unfiltered"]:
        do_reports(numerical_filter, exp_1, exp_2, label_width)
    else:
        do_reports(lambda x:x, exp_1, exp_2, label_width)

    