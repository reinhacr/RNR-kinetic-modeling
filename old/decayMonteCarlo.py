import argparse
import numpy as np
import pandas as pd
from run_gillespie_fx import runToSteadyState, exponentialFitGillespie, runGillespie
from gillespieCore_fx import parseGillespie
from myTools.misc import readSimpleInput
from tqdm import tqdm


def decayMonteCarlo(*, param, monitoredIdx, max_steps, averagingSteps, tol,
                    minRSquared,
                    nSteps, highTemp, lowTemp, coolingSchedule, tauWanted, rateStep,
                    logfile, **kwargs):

    if lowTemp < 1e-10:
        lowTemp = 1e-20
    # param.counts is changed in-place for speed
    # need to reset it every time
    originalCounts = param.counts.copy()
    time, counts = runToSteadyState(param, monitoredIdx, max_steps, averagingSteps, tol)
    param.counts = originalCounts

    tau = exponentialFitGillespie(time, counts, minRSquared=minRSquared)
    costFunction = np.abs(tau - tauWanted)

    if pd.isna(costFunction):
        costFunction = 1e10

    if coolingSchedule == "exponential":
        decay = (1 / nSteps) * np.log(highTemp / lowTemp)
        temperatures = [highTemp * np.exp(-decay * i) for i in range(nSteps)]
    else:
        # linear cooling schedule
        slope = (highTemp - lowTemp) / nSteps
        temperatures = [highTemp - slope * i for i in range(nSteps)]

    energy = costFunction
    allEnergies = np.zeros(nSteps)
    nAccepted = 0

    log = open(logfile, "w+")
    for i in tqdm(range(len(temperatures))):

        oldRateConstants = param.rate_constants.copy()
        param.rate_constants = newRateConstants(param, rateStep)  # implement constraints here

        # param.counts is changed in-place during the simulation (for speed)
        # need to reset it every time
        originalCounts = param.counts.copy()
        time, counts = runToSteadyState(param, monitoredIdx, max_steps, averagingSteps, tol)
        param.counts = originalCounts

        tau = exponentialFitGillespie(time, counts, minRSquared=minRSquared)

        if pd.isna(tau) is True or np.isinf(tau) is True:
            tau = 1e10  # so it will "never" be accepted

        costFunctionNew = np.abs(tau - tauWanted)
        if pd.isna(costFunctionNew):
            costFunctionNew = 1e10

        if pd.isna(costFunction):
            costFunction = 1e10

        beta = 1 / temperatures[i]
        if costFunctionNew < costFunction:
            energy = costFunctionNew
            costFunction = costFunctionNew
            nAccepted += 1
        elif np.random.rand() < np.exp(-beta * (costFunctionNew - costFunction)):
            energy = costFunctionNew
            costFunction = costFunctionNew
            nAccepted += 1
        else:
            param.rate_constants = oldRateConstants

        allEnergies[i] = energy

        # print("step:", i, "of", nSteps, " temperature:", temperatures[i], " energy:", energy)
        log.write(f"step: {i} of {nSteps}, temperature: {temperatures[i]}, energy: {energy}\n")
        log.flush()
        print(f"step: {i} of {nSteps}, temperature: {temperatures[i]}, energy: {energy}\n", flush=True)
    log.close()
    return param.rate_constants, allEnergies, nAccepted / nSteps, tau

def newRateConstants(param, rateStep):

    ks = param.rate_constants.copy()

    goAhead = False
    conditions = []
    conditions.append(False)
    conditions.append(False)
    conditions.append(False)
    conditions.append(False)
    # conditions.append(False)

    while goAhead is False:
        # print("trying to satisfy conditions...")
        newks = ks.copy()  # trial ks
        rxnIdx = np.random.randint(param.n_reactions - 2)
        # shift = np.random.choice([-1, 1]) * rateStep
        shift = np.random.randn() * rateStep
        newks[rxnIdx] += shift

        # if rxnIdx % 2 == 0:
        #     newks[rxnIdx + 1] -= shift
        # elif rxnIdx % 2 == 1:
        #     newks[rxnIdx - 1] -= shift

        for i in range(len(newks)):
            if newks[i] < 0.0:
                newks[i] = 0.0

        # Y731 <-> Y730 roughly isoergic
        if 0.8 < newks[4] / newks[5] and newks[4] / newks[5] < 1.2:
            conditions[0] = True

        # Y122 -> Y356 is approximately zero for photoRNR
        if newks[1] < 1e-6:
            conditions[1] = True

        # Y731 <-> Y730 > Y730 <-> C439
        if newks[4] + newks[5] > newks[6] + newks[7]:
            conditions[2] = True
        # goAhead = True

        if 31250 * 0.95 < np.sum(newks[0:8]) < 31250 * 1.05:
            conditions[3] = True

        # # forward rates < backward rate constants
        # if ks[1] + ks[3] + ks[5] + ks[7] < ks[0] + ks[2] + ks[4] + ks[6]:
        #     conditions[4] = True

        if all(conditions) is True:
            ks = newks
            goAhead = True

    # print(list(ks))
    return ks


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    mainArgs = ap.parse_args()
    args = readSimpleInput(mainArgs.infile)
    args.param = parseGillespie(args.infile)

    print("running simulated annealing...")
    rate_constants, allEnergies, fractionAccepted, tau = decayMonteCarlo(**vars(args))

    args.param.rate_constants = rate_constants

    with open(args.out + ".csv", "w+") as csv:
        for rxn, k in zip(args.param.reaction, args.param.rate_constants):
            csv.write(rxn)
            csv.write(",")
            csv.write(str(k))
            csv.write(",\n")

    print("running simulation with new rate constants...")
    runGillespie(args.final_simulation_time, args.param, args.out + "_final_rxn.csv")

    np.savetxt(args.out + "_energies.csv", allEnergies.T, delimiter=",")

    with open(args.out + "_fraction_accepted.csv", "w+") as csv:
        csv.write(str(fractionAccepted))
        csv.write(" accepted,")

    import matplotlib.pyplot as plt
    plt.plot(allEnergies)
    plt.title("energy")
    plt.show()
