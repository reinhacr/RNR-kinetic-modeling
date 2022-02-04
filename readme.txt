The csv file is provided as test1_1.csv (I usually generate an array of these, hence the numbering). The reduced species have been encorporated in, so you now have 2 reactants and four terms in the stoich matrix for each case, and i added two "Dummy" species to keep the y731 flipping reaction with the same number of reactants as everything else. 
The following is what i use to run this in the dsq job list
module load Langs/Python/3.6-anaconda; source activate plot_env; python decayMonteCarlo_full.py test1_1.in > test1_1.txt
Here, I have my numba and other installations in a conda environment. You may have something else or choose to run this on your computer
