import matplotlib

font = {'family' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt

import csv
import pandas as pd

df = pd.read_csv('result.csv')

plt.plot(df["num_cons"], df["time"], "kx--")

plt.xlabel("Execution time (seconds)")
plt.ylabel("Number of services")
plt.title("Execution time for varying number of services")

plt.savefig("figureTimePRES.png", dpi=600)