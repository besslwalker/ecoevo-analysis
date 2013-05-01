#!/bin/tcsh

# create the file anew
rm lineage_test_raw.dat
touch lineage_test_raw.dat

# for each run, list the line count, file name, and last line of the lineage_domEQU.dat file
foreach i (limited-81_7???)
  wc -l $i/data/lineage_domEQU.dat >> lineage_test_raw.dat
  tail -1 $i/data/lineage_domEQU.dat >> lineage_test_raw.dat
end

# call a python script to parse all that information into a human-skimmable format

python lineage_test_consolidator.py
