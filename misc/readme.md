# This directory contains additional scripts and advice for parsing output of the kinetic model
Contents:

template.csv:

A template csv file where placeholder strings for k1, k2, k3, and k4 are placed. csv_prep.sh will then calculate the effective rate constant expected using the relationship with keq_conf (flipped v. stacked equillibrium) and update your csv file with the correct rate constant. This step is handled by a bash script so that the code stays very general and if users want to apply a scalar correction to rate constants they can do so using the provided bash script. 

csv_prep.sh:

A bash script to calculate k1eff,k2eff,k3eff,and k4eff and replace the placeholder strings in the template.csv file. Python or any other language would also work fine for this task as its a simple script. 

rnr_rate_constants_literature_comparison.xlsx:

Not peer reviewed or included in any published work, but collection of rate constants from literature with sources, and the interpretation of what I believe they are measuring gathered in the initiation of this work. I am providing it here in case it is useful to the community at some point. Please email clorice.reinhardt@aya.yale.edu if you have any comments on this or want to submit an addition or correction. I am open to expanding this to more recent studies or RNRs from other species as well if someone wants to help maintain or grow it. 

General advice:

The code as deposited is very useful. If you want to obtain plots for different radical amounts over time, extend the time for which the kinetic model is solved beyond the shorter time periods used to compare to experimental data and add additional monitored indexes. This is generally not recommended to do during Monte Carlo/Simulated Annealing runs and initial exploration of rate constant space as it will significantly increase the time required for a calculation. First, find a desired set of rate constants, inspect the exponential fit and Y356 profile (if you are sticking true to our original photoRNR model), then solved for other radical amounts over time while not making use of the Monte Carlo/Simulated Annealing framework. 

If extending amount of time that rate matrix is solved, please adjust nSteps so resolution can be tuned and you can still get enough time points out to suite your purpose
