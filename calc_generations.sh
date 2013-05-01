#!/bin/tcsh

# Uses tail -1 to create a file with the last lines of each time.dat in the 
# experiment.
#
# BLW
# 9-27-11

# create the file anew
rm time_tails.dat
touch time_tails.dat

foreach i ( unlimited-base_1* limited-base_2* severe-base_3* unlimited-21_4* unlimited-81_6* unlimited-3task_8* limited-3task_9* severe-3task_10* unlimited-noEQU_11* limited-noEQU_12* severe-noEQU_13* unlimited-noloss_14* limited-noloss_15* )
	gunzip $i/data/time.dat.gz
	tail -1 $i/data/time.dat >> time_tails.dat
	gzip --best $i/data/time.dat
end

# call python script to do percentile/median reporting
python calc_generations.py