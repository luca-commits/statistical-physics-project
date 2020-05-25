import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# from numpy import genfromtext

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def main():
    data_fullrun = np.genfromtxt('energies_fullrun.out')
    data_nobonds = np.genfromtxt('energies_nobonds.out')
    data_nobonds_noangles = np.genfromtxt('energies_nobonds_noangles.out')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5), dpi=300)
    fig.tight_layout(pad=5)

    ax1.set_yscale('log')
    ax1.minorticks_on()

    ax2.set_yscale('log')
    ax2.minorticks_on()

    ax2.plot(data_fullrun[:, 1], data_fullrun[:, 2], label='full run')
    ax2.legend()

    ax1.plot(data_nobonds[:, 1], data_nobonds[:, 4], label='bond terms off')
    ax1.plot(data_nobonds_noangles[:, 1], data_nobonds_noangles[:, 4], label='bond and angle terms off')
    ax1.legend()

    plt.savefig('plot.png')
    plt.show()

if __name__ == "__main__":
  main()