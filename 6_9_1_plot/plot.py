import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def plotCO(ax):
  # For C-O bond, values taken from Cornell (1995) paper (p. 5194)and Morse (1929) paper (p. 62)
  D_e = 5.65e7 # 91600 # wave no = 1 / cm (online)
  r_0 = 1.229 # Angström (Morse)
  a = 2.29 # 1 / Angström (Morse)
  k_r = 0.06 # 12.53e-3 # 1 / (Angström^3) // // 1900 N / m (online)

  r = np.arange(0, 4, 0.01)

  V_harmonic = 0.5 * k_r * (r - r_0)**2 * 10e10 # wave no = 1 / cm
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2 # wave no = 1 / cm

  ax.plot(r, V_morse, label='Morse-Potential')
  ax.plot(r, V_harmonic, label='harmonisches Potential')
  ax.axhline(D_e, 0, 1, color='black', linestyle='--', linewidth=1, label='Dissoziierungsenergie')

  ax.set_ylim(0, 10e7)
  ax.set_xlim(0, 4)

  ax.minorticks_on()

  # ax.set(title="$ CO $ bond")
  ax.set(xlabel="Abstand der Atome (in Angström)")
  ax.set(ylabel="Energie (in $ 1 / m $)")
  ax.legend()

def plotN2(ax):
  # For N_2 bond, values taken from Cornell (1995) paper (p. 5194)and Morse (1929) paper (p. 62)
  D_e = 4.95e7 # wave no = 1 / cm // 945 kJ / mol (online)
  r_0 = 1.09 # Angström (Morse)
  a = 2.56 # 1 / Angström (Morse)
  k_r = 0.0723 # 1 / (Angström^3) // k_r = 2287 N/m (online)

  r = np.arange(0, 4, 0.01)

  V_harmonic = 0.5 * k_r * (r - r_0)**2 * 10e10 # wave no = 1 / m
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2 # wave no = 1 / m

  ax.plot(r, V_morse, label='Morse potential')
  ax.plot(r, V_harmonic, label='harmonic potential')
  ax.axhline(D_e, 0, 1, color='black', linestyle='--', linewidth=1, label='dissociation energy')

  ax.set_ylim(0, 10e7)
  ax.set_xlim(0, 4)

  ax.minorticks_on()

  ax.set(title="$ N_2 $ bond")
  ax.set(xlabel="nuclear seperation (in Angström)")
  ax.set(ylabel="energy (in $ 1 / m $)")
  ax.legend()

def main():
  # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8,7.2), dpi=300)
  fig, ax1 = plt.subplots(1, 1, figsize=(5,5), dpi=600)
  fig.tight_layout(pad=2.5)

  plotCO(ax1)
  # plotN2(ax2)

  plt.suptitle("Vergleich von Näherungen für das Potential für ein $ CO $-Molekül", fontweight='bold')

  plt.savefig("plot.png")

if __name__ == "__main__":
  main()
