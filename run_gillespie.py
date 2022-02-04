import pandas as pd
import matplotlib.pyplot as plt
from run_gillespie_fx import runGillespie
from gillespieCore_fx import parseGillespie

infile = "test_unimolecular.csv"
outfile = "test_gillespie.csv"
time = 0.01  # seconds
saveFreq = 10  # steps

param = parseGillespie(infile)
runGillespie(time, param, outfile, saveFreq)

# df = pd.read_csv(outfile)
# plt.plot(df["time"], df["Y356_dot"])
# plt.xlabel("time (s)")
# plt.ylabel("Y356• level")
# plt.show()
