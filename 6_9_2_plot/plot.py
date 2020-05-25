import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# from numpy import genfromtext

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def main():
    data_fullrun = np.genfromtxt('energies_fullrun.out')
    data_noangles = np.genfromtxt('energies_noangles.out')

    fig, ax = plt.subplots(1, 1, figsize=(5,5), dpi=300)
    fig.tight_layout(pad=2)

    ax.minorticks_on()
    ax.set(xlabel="Zeit (in $ ps $)")
    ax.set(ylabel="potentielle Energie (in $ kJ / mol $)")

    ax.plot(data_fullrun[:, 1], data_fullrun[:, 4], label='LJ, bond, angle, dihedral terms')
    ax.plot(data_noangles[:, 1], data_noangles[:, 4], label='LJ, bond, dihedral terms')
    ax.legend()

    plt.savefig('plot.png')

if __name__ == "__main__":
  main()
