#!/bin/bash
# Put Keq here
# I construct a template input file, and the values are computed with bash and bc (because of floating point operations). 7 decimal places are kept. This could be easily worked into a python workflow, but works well enough
# Note in input files for KM, E can be used in scientific notation, but for bash it needs to be written out explicitly, i.e instead of 5.2e9, 5.2x10^9 is needed
# echo statements used liberally for sanity checks
Keq=$(bc <<<"scale=7;1")
k1=$(bc <<<"scale=7;86000")
k2=$(bc <<<"scale=7;8600")
k3=$(bc <<<"scale=7;5.2*10^9")
k4=$(bc <<<"scale=7;6100")
echo $k1 
echo $k2
echo $k3
echo $k4
echo $Keq
keq1=$(bc <<<"scale=7;1/($Keq + 1)")
echo $keq1
keq2=$(bc <<<"scale=7;$Keq/($Keq + 1)")
echo $keq2
k1_eff=$(bc <<<"scale=7;$k1*$keq1")
k2_eff=$(bc <<<"scale=7;$k2*$keq1")
k3_eff=$(bc <<<"scale=7;$k3*$keq2")
k4_eff=$(bc <<<"scale=7;$k4*$keq2")
echo $k1_eff
echo $k2_eff
echo $k3_eff
echo $k4_eff
# Input effective rate constants into spreadsheet for km
sed -e 's/k1/'"$k1_eff"'/ ; s/k2/'"$k2_eff"'/ ; s/k3/'"$k3_eff"'/; s/k4/'"$k4_eff"'/' template.csv > mod_keq.csv
#done
