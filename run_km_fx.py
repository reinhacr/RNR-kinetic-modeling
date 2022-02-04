import os
import sys
import numpy as np
import pandas as pd
import scipy.optimize as optimize
from scipy.linalg import expm
from km_core_fx import gillespieStep
import matplotlib.pyplot as plt


def runGillespie(time, param, outfile, saveFreq=1):

    if os.path.exists(outfile):
        os.remove(outfile)

    print("setting initial conditions...")
    outfile = open(outfile, "ab")
    t = 0.0
    step = 0

    print("starting simulation...")
    np.savetxt(outfile, [], header="time,rxnIdx," + ",".join(param.speciesNames), comments="")

    # IMPORTANT CHECK TO PREVENT SEGFAULTS AND FAILURE OF THE ALGORITHM
    rates = np.zeros(param.n_reactions)

    for i in range(param.n_reactions):
        rates[i] = param.rate_constants[i]
        for j in range(param.n_reactants[i]):
            rates[i] *= param.counts[param.species[i][j]]
    if np.sum(rates) == 0.0:
        sys.exit("ERROR: initial total rate is zero, system cannot move from here")

    while t < time:
        t, reaction_idx, param.counts = gillespieStep(t, param.n_reactions,
                                                      param.counts, param.rate_constants,
                                                      param.species, param.stoich,
                                                      param.n_species_involved, param.n_reactants)
        step += 1
        if (step % saveFreq) == 0:
            output = np.concatenate(([t, reaction_idx], param.counts))
            np.savetxt(outfile, output.reshape(1, output.shape[0]), delimiter=",", fmt="%g")

        if step % 1000 == 0:
            print(f"{100 * t / time:0.1f}% complete", end="\r")
    outfile.close()
    print("simulation done")


def runGillespieInMemory(time, param):
    """this bad boy is very slow because it keeps having to allocate memory for the growing data"""
    t = 0.0
    step = 0

    output = pd.DataFrame(columns=["time", "reaction_idx"] + list(param.speciesNames))
    while t < time:
        t, reaction_idx, param.counts = gillespieStep(t, param.n_reactions,
                                                      param.counts, param.rate_constants,
                                                      param.species, param.stoich,
                                                      param.n_species_involved, param.n_reactants)
        output.loc[step] = [time, reaction_idx] + list(param.counts)
        step += 1

    return output


def loadGillespieSimulation(infile):
    data = pd.read_csv(infile, engine="c")
    return data


def runToSteadyState(param, monitored_idx=0, max_steps=100000, averaging_steps=10000, tol=1e-1):
    """
    run a Gillespie simulation until the concentration of the indicated species stabilizes
    to get new counts

    Arguments
    ---------
    max_steps: int
    averaging_steps: int
        how many steps to sample to determine stabilizaion
    tol: float
        stabilization tolerance, default=1e-2
    monitored_idx: int
        which species to follow to determine stabilization
    """

    t = 0.0
    monitored_counts = np.zeros(max_steps) * np.nan
    times = np.zeros(max_steps) * np.nan

    step = 0
    oldStep = 0

    counts = param.initial_counts.copy()
    for i in range(max_steps):
        # always start with the same counts
        t, reaction_idx, counts = gillespieStep(t, param.n_reactions,
                                                counts, param.rate_constants,
                                                param.species, param.stoich,
                                                param.n_species_involved, param.n_reactants)

        monitored_counts[step] = counts[monitored_idx]
        times[step] = t
        step += 1
        if i - oldStep == averaging_steps * 2:
            halfWayStep = i - averaging_steps

            firstMean = np.mean(monitored_counts[oldStep:halfWayStep])
            secondMean = np.mean(monitored_counts[halfWayStep:i])

            if np.abs(secondMean - firstMean) < tol and np.abs(secondMean - firstMean) != 0.0:
                break
            else:
                oldStep = i
    else:
        # print(f"warning: hit max iterations without converging monitored_idx {monitored_idx}!\n")
        pass

    # plt.plot(times, monitored_counts)
    # plt.show()
    return times[~np.isnan(times)], monitored_counts[~np.isnan(monitored_counts)]


def exponentialFitTrajectory(times, counts, minRSquared=0.95):
    # data = pd.read_csv(infile, engine="c")
    # time = data["time"].values
    # counts = data[speciesToFit]

    # fit to 3 microseconds onward, like in experiment
    mask = np.logical_and(times > 3.0E-6, times < 90.0E-6)

    # times = times[mask]
    # counts = counts[mask]

    # import matplotlib.pyplot as plt

    # take the first half of the decay to remove the flat part (can cause fit to fail)
    # time = time[0:len(time) // 2]
    # counts = counts[0:len(counts) // 2]

    plt.plot(times, counts)
    plt.show()

    def monoExp(x, m, k, b):
        return m * np.exp(-k * x) + b

    p0 = (0.0, 0.0, 0.0)  # guess

    try:
        params, cv = optimize.curve_fit(monoExp, times, counts, p0)
        m, k, b = params

        tau = 1 / k
        # squaredDiffs = np.square(counts - monoExp(time, m, k, b))
        # squaredDiffsFromMean = np.square(counts - np.mean(counts))
        # rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
        # print(rSquared)
        # if rSquared < minRSquared:
        #     raise RuntimeError(f"exponential fit failed: rSquared = {rSquared} is too small")
        # plt.plot(time, monoExp(time, m, k, b), color='r',label='fit: m=%5.3f, k=%5.3f, b=%5.3f' % tuple(params))
        # plt.plot(time, counts, color='b', label='Y356 counts')
        # plt.ylabel("counts")
        # plt.xlabel("time")
        # plt.legend()
        # plt.show()
        # plt.plot(time, monoExp(time, m, k, b))
        # plt.show()

        # fit = np.polyfit(time, np.log(counts), 1)
        # plt.plot(time, np.log(counts))
        # print(fit)
        # print(tau)
        # tau = 1 / -fit[0]
    except RuntimeError:
        tau = np.nan

    return tau


def solveRateMatrix(param, monitoredIdx, solveMatrixTime):
    # param includes rate_constants
    # build rate matrix

    nSpecies = len(param.initial_counts)
    nRxns = len(param.rate_constants)
    simpleRateMatrix = np.zeros((nSpecies, nSpecies))

    for i in range(nSpecies):
        for j in range(nSpecies):
            for k in range(nRxns):
                if i == param.species[k][0] and j == param.species[k][1] and i != j:
                    simpleRateMatrix[i, j] = param.rate_constants[k]

    odeMatrix = np.zeros((nSpecies, nSpecies))
    for i in range(nSpecies):
        for j in range(nSpecies):
            if i == j:
                for k in range(nSpecies):
                    if i != k:
                        odeMatrix[i, j] += simpleRateMatrix[i, k]
            else:
                odeMatrix[i, j] = -simpleRateMatrix[j, i]

    nSteps = 400

    times = np.array([(i / nSteps) * solveMatrixTime for i in range(nSteps)])
    counts = np.array([(expm(-odeMatrix * t) @ param.initial_counts)[monitoredIdx] for t in times])

    #print(times)
    #print(counts)
    # plt.plot(times, counts)
    # plt.show()

    return times, counts


# if __name__ == "__main__":
#     infile = "test_unimolecular.csv"
#     outfile = "test.csv"
#     param = parseGillespie(infile)

#     # start = time()
#     # runGillespie(1e-3, param, outfile, saveFreq=1)
#     # print(time() - start)

#     # start = time()
#     # df = loadGillespieSimulation(outfile)
#     # print(time() - start)

#     time, counts = runToSteadyState(param, outfile, monitored_idx=1, max_steps=100000000, averaging_steps=10000, tol=1e-1)
#     tau = exponentialFitGillespie(time, counts)
