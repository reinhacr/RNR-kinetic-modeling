import numpy as np
import numba
from types import SimpleNamespace
import pandas as pd


def parseInput(infile):

    df = pd.read_csv(infile)
    headers = df.columns
    named_headers = [header for header in headers if not header.startswith("Unnamed")]

    headerStr = " ".join(named_headers)
    headerShouldBe = " ".join(["species_idx", "species_names", "initial_counts", "counts_max",
                               "reaction", "reaction_idx", "rate_constants", "species", "stoich", "n_reactants"])

    if not headerStr.startswith(headerShouldBe):
        raise ValueError("headers must be "
                         "species_idx, species_names, initial_counts, counts_max, reaction, reaction_idx, rate_constants, species, stoich, n_reactants")

    speciesInfo = df[["species_idx", "species_names", "initial_counts", "counts_max"]].dropna()
    rxnInfo = df.drop(["species_idx", "species_names", "initial_counts", "counts_max"], axis=1).dropna()

    param = SimpleNamespace()
    param.speciesNames = speciesInfo.species_names.values

    param.initial_counts = np.array(speciesInfo.initial_counts.values).astype(np.int)
    param.counts_max = np.array(speciesInfo.counts_max.values).astype(np.int)
    param.rate_constants = np.array(rxnInfo.rate_constants.values).astype(np.float64)
    param.n_species = len(speciesInfo.species_names)
    param.n_reactants = np.array(rxnInfo.n_reactants).astype(np.int)
    param.n_reactions = len(param.rate_constants)

    headers = rxnInfo.columns
    named_headers = [header for header in headers if not header.startswith("Unnamed")]

    speciesColumns = []
    col = 0
    for i in range(len(headers)):
        col = i
        if headers[col] == "species":
            break

    for i in range(col, len(headers)):
        col = i
        if headers[col] == "stoich":
            break
        else:
            speciesColumns.append(rxnInfo[headers[col]].values)

    stoichColumns = []
    for i in range(col, len(headers)):
        col = i
        if headers[col] == "n_reactants":
            break
        else:
            stoichColumns.append(rxnInfo[headers[col]].values)

    param.species = np.array(speciesColumns).T.astype(np.int)
    param.stoich = np.array(stoichColumns).T.astype(np.int)

    if param.species.shape != param.stoich.shape:
        raise ValueError("species and stoich must have same shape")

    param.n_species_involved = np.array([len(row[~np.isnan(row)]) for row in param.species]).astype(np.int)
    param.species[np.isnan(param.species)] = -1
    param.stoich[np.isnan(param.species)] = -1

    # make all arrays memory-contiguous for speed
    dictParam = vars(param)
    for key, value in dictParam.items():
        if type(dictParam[key]) is np.ndarray:
            dictParam[key] = np.ascontiguousarray(dictParam[key])

    param.reaction = df.reaction.values

    return param


@numba.jit(nopython=True)
def gillespieStep(t, n_reactions, counts, rate_constants, species, stoich, n_species_involved, n_reactants):

    rates = np.zeros(n_reactions)

    # for bimolecular reactions
    # goAhead = False
    # while goAhead is False:
    #     for i in range(n_reactions):
    #         rates[i] = rate_constants[i]
    #         for j in range(n_reactants[i]):
    #             rates[i] *= counts[species[i, j]]   # for unimolecular reactions, concentration in molar, rates in molar/s

    #     totalRate = np.sum(rates)
    #     if totalRate == 0.0:
    #         print("gillespieStep(): system reached a stationary state. Perturbing to get it moving again...")
    #         for i in range(len(counts)):
    #             if counts[i] == 0:
    #                 counts[i] += 1
    #                 if i % 2 == 0:
    #                     counts[i + 1] -= 1
    #                 else:
    #                     counts[i - 1] -= 1
    #     else:
    #         goAhead = True

    for i in range(n_reactions):
        rates[i] = rate_constants[i]
        for j in range(n_reactants[i]):
            rates[i] *= counts[species[i, j]]   # for unimolecular reactions, concentration in molar, rates in molar/s
    totalRate = np.sum(rates)

    if totalRate == 0.0:
        return t + 0, -1, counts

    rates /= totalRate
    # make dartboard
    dartboard = np.zeros(n_reactions + 1)
    for i in range(n_reactions):
        for j in range(i + 1):
            dartboard[i + 1] += rates[j]

    rxnIdx = 0
    rand01 = np.random.rand()
    for i in range(n_reactions):
        if dartboard[i] < rand01 and rand01 < dartboard[i + 1]:
            rxnIdx = i
            break

    for i in range(n_species_involved[rxnIdx]):
        counts[species[rxnIdx, i]] += stoich[rxnIdx, i]

    rand01 = np.random.rand()
    randTime = (1 / totalRate) * np.log(1 / rand01)

    t += randTime

    return t, rxnIdx, counts
