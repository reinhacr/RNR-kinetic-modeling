import argparse
import numpy as np
import pandas as pd
from run_km_fx import runToSteadyState, exponentialFitTrajectory, runGillespie, solveRateMatrix
from km_core_fx import parseInput
from myTools.misc import readSimpleInput
from tqdm import tqdm


def decayMonteCarlo(*, param, monitoredIdx, max_steps, averagingSteps, tol,
                    minRSquared, solveMatrixTime=None,
                    nSteps, highTemp, lowTemp, coolingSchedule, plateauTemp=None, tauWanted, rateStep,
                    logfile, **kwargs):

    if lowTemp < 1e-10:
        lowTemp = 1e-20

    # originalCounts = param.counts.copy()  # runToSteadyState changes counts, see below

    # --- this part can either be gillespie or rate matrix exponentiation ---

    times, counts = runToSteadyState(param, monitoredIdx, max_steps, averagingSteps, tol=10**-10)
    np.savetxt("gillespieRun.csv", np.array([times, counts]).T, delimiter=",")
    #print(param.initial_counts)
    times2, counts2 = solveRateMatrix(param, monitoredIdx, solveMatrixTime)
    np.savetxt("solveRateMatrix.csv", np.array([times2, counts2]).T, delimiter=",")
    # print(param)
    # print(times)
    # print(counts)
    # exit()
    #exit()
    import matplotlib.pyplot as plt
    plt.plot(times, counts, label="Gillespie")
    plt.plot(times2, counts2, label="Matrix Exp.")
    plt.xlim([0, 4e-5])
    plt.xlabel("Time (s)")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()
    exit()
    # -----------------------------------------------------------------------
    #param.counts = originalCounts  # runToSteadyState changes counts, see above

    tau = exponentialFitTrajectory(times, counts, minRSquared=minRSquared)
    print(tau)
    exit()



    costFunction = np.abs(tau - tauWanted)

    if pd.isna(costFunction):
        costFunction = 1e10

    if coolingSchedule == "exponential":
        decay = (1 / nSteps) * np.log(highTemp / lowTemp)
        temperatures = [highTemp * np.exp(-decay * i) for i in range(nSteps)]
    elif coolingSchedule == "linear":
        # linear cooling schedule
        slope = (highTemp - lowTemp) / nSteps
        temperatures = [highTemp - slope * i for i in range(nSteps)]
    elif coolingSchedule == "linear plateau":
        slope = (highTemp - lowTemp) / (nSteps // 2)
        temperaturesSlope = [highTemp - slope * i for i in range(nSteps // 2)]
        temperaturePlateau = [plateauTemp for i in range(nSteps // 2)]
        temperatures = temperaturesSlope + temperaturePlateau

    energy = costFunction
    allEnergies = np.zeros(nSteps)
    nAccepted = 0

    log = open(logfile, "w+")
    for i in tqdm(range(len(temperatures))):

        oldRateConstants = param.rate_constants.copy()
        param.rate_constants = newRateConstants(param, rateStep)  # implement constraints here

        originalCounts = counts.copy()  # runToSteadyState changes counts, see below

        # --- this part can either be gillespie or rate matrix exponentiation ---
        time, counts = runToSteadyState(param, monitoredIdx, max_steps, averagingSteps, tol)
        # times, counts = solveRateMatrix(param, monitoredIdx, solveMatrixTime)
        # ------------------------------------------------------------------------
        param.counts = counts  # runToSteadyState changes counts, see above

        tau = exponentialFitTrajectory(times, counts, minRSquared=minRSquared)

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
        tau_print = tau * 1e6
        # print("step:", i, "of", nSteps, " temperature:", temperatures[i], " energy:", energy)
        log.write(f"step: {i} of {nSteps}, temperature: {temperatures[i]}, energy: {energy}\n")
        log.flush()
        print((param.rate_constants.tolist()), tau_print)
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
    conditions.append(False)
    conditions.append(False)

    while goAhead is False:
     #   print("trying to satisfy conditions...")
        newks = ks.copy()  # trial ks
        rxnIdx = np.random.randint(param.n_reactions - 5)
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

        # Y731 <-> Y730 is 0 in y730f
#        if newks[4] < 1e-6:
#            conditions[0] = True

        # Y356 -> Y122 is approximately zero for photoRNR
        if newks[0] < 1e-6:
            conditions[0] = True

        # Y122 -> Y356 is approximately zero for photoRNR
        if newks[1] < 1e-6:
            conditions[1] = True

        # Y356 -> Re is approximately zero for photoRNR
        if newks[9] < 1e-6:
            conditions[2] = True

        # Sink -> Y356 is approximately zero for photoRNR
        if newks[11] < 1e-6:
            conditions[3] = True

        # Radical decay rate to sink should be less than or equal to 55000 but not 0
        if 1 * 0.98 < newks[6] < 55000 * 1.02:
            conditions[4] = True

        # C439 <-> Y730 is 0 in y730f
#        if newks[7] < 1e-6:
#            conditions[5] = True

        # Y731 <-> Y730 is 0 in y730f
#        if newks[5] < 1e-6:
#            conditions[7] = True

        # C439 <-> Y730 is 0 in y730f
#        if newks[8] < 1e-6:
#            conditions[6] = True

        # Enforce K_eq around 3-6
#        if 3.2 < (newks[3] / newks[2]) < 6.2:
#            conditions[9] = True

        # RT to Y356 from Y731 and from Y356 to 731 need to be greater then enzyme turnover in wild-type RNR
        if newks[2] and newks[3] > 10:
            conditions[5] = True

        if all(conditions) is True:
            ks = newks
            goAhead = True
#    print(list(ks))
    return ks


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    mainArgs = ap.parse_args()
    args = readSimpleInput(mainArgs.infile)
    args.param = parseInput(args.infile)

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
    from matplotlib import rc
    rc("font", **{"family": 'sans-serif', "sans-serif": ['Arial'], "size": 14})
    plt.plot(allEnergies)
    plt.ylabel("Simulated Annealing Energy", fontsize=14)
    plt.xlabel("Steps", fontsize=14)
    #   plt.title("energy")
    plt.tight_layout()
    plt.show()

    df = pd.read_csv("tauWanted_decay_730f_final_rxn.csv")
    plt.plot(df["time"], df["Y356_dot"])
    plt.xlabel("time (s)")
    plt.ylabel("Y356• level")
    plt.show()
