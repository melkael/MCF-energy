import matplotlib

font = {'family' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt

import csv
import pandas as pd

df = pd.read_csv('result.csv')

plt.plot(df["num_cons"], df["obj"], "kx--")

plt.xlabel("Energy consumption (Watts)")
plt.ylabel("Number of services")
plt.title("Energy consumption for varying number of services")

plt.savefig("figureEnergyPRES.png", dpi=600)