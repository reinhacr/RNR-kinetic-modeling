# RNR-kinetic-modeling

In the current form, this code allows you to sample sets of rate constants for a process such as radical transfer in RNR keeping in mind an arbitrary set of constraints and a cost function consisting of the absolute value of the difference between some desired decay rate and the rate of decay produced by the current set of rate constants. This code can easily be modified to explore rate constants for something other than RNR as long as

1) you can define a suitable cost function like the rate of decay of a species to guide your sampling
2) the process you are studying involves first-order kinetics
3) you understand how the input files work

The code is structured as follows:
the main module is called **decayMonteCarlo_full.py**
  this controls all the other code and contains the function that runs the Monte Carlo sampling of rate constants. This is where you can define constraints on the rate constants (see, for example, the constaints set near line 90) and the cost function (line 24).

this module depends on the modules **run_km_fx.py** and **km_core_fx.py**, as well as the directory **myTools**, which contains the function that reads the main (text) input file. **run_km_fx.py** contains functions to run Gillespie simulations of dynamics, perform exponential fitting of radical decay to evaluate the cost function, run Gillespie simulations to steady state, and to solve the rate matrix (matrix exponential solution to a system of ODEs). **km_core_fx.py** contains functions to parse the csv input file, and take a Gillespie step (numba-compiled for speed).

The code requires **two** input files, one "text" input file and one csv input file. The text input file contains a path to the csv file and various arguments. The csv input file contains the full description of the chemical reactions being studied. The arguments in the text input file are:

infile: [your csv file]

out: [outfile]

logfile: [logfile]

### gillespie and fitting parameters
monitoredIdx: [index of chemical species to monitor (from csv)]

max_steps: [max number of steps to take when seeking steady state]

averagingSteps: [number of steps to use at a time when determining steady state]

tol: [tolerance to use when determining steady state]

minRSquared: [minimum r^2 to accept when fitting a exponential]

### monte carlo parameters
nSteps: [number of Monte Carlo steps]

highTemp: [initial "temperature" for simulated annealing]

lowTemp: [final "temperature" for simulated annealing, can be set to same as highTemp for steady "temperature"]

coolingSchedule: [linear or exponential]

tauWanted: [desired decay rate]

rateStep: [how much to vary rate constants between Monte Carlo steps]

final_simulation_time: [length of Gillespie simulation]
