import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def plotCO(ax):
  # For C-O bond, values taken from Cornell (1995) paper (p. 5194)and Morse (1929) paper (p. 62)
  D_e = 91600 # wave no = 1 / cm
  r_0 = 1.229 # Angström
  a = 2.29 # 1 / Angström
  k_r = 12.53e-3 # 1 / (Angström^3)
  
  r = np.arange(0, 4, 0.01)
  
  V_harmonic = 0.5 * k_r * (r - r_0)**2 * 10e8 # wave no = 1 / cm
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2 # wave no = 1 / cm
  
  ax.plot(r, V_morse, label='Morse potential')
  ax.plot(r, V_harmonic, label='harmonic potential')
  ax.axhline(D_e, 0, 1, color='black', linestyle='--', linewidth=1, label='dissociation energy') 
  
  ax.set_ylim(0, 150000)
  ax.set_xlim(0, 4)
  
  ax.minorticks_on()
  
  ax.set(title="$ CO $ bond")
  ax.set(xlabel="nuclear seperation (in Angström)")
  ax.set(ylabel="energy (in wave numbers)")
  ax.legend()

def plotN2(ax):
  # For H_2 bond, values taken from Cornell (1995) paper (p. 5194)and Morse (1929) paper (p. 62)
  D_e = 40100 # wave no = 1 / cm
  r_0 = 0.76 # Angström
  a = 2.56 # 1 / Angström
  k_r = 16.131e-3 # 1 / (Angström^3) // k_r = 510 N/m
  
  r = np.arange(0, 4, 0.01)
  
  V_harmonic = 0.5 * k_r * (r - r_0)**2 * 10e8 # wave no = 1 / cm
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2 # wave no = 1 / cm
  
  ax.plot(r, V_morse, label='Morse potential')
  ax.plot(r, V_harmonic, label='harmonic potential')
  ax.axhline(D_e, 0, 1, color='black', linestyle='--', linewidth=1, label='dissociation energy') 
  
  ax.set_ylim(0, 80000)
  ax.set_xlim(0, 4)
  
  ax.minorticks_on()
  
  ax.set(title="$ N_2 $ bond")
  ax.set(xlabel="nuclear seperation (in Angström)")
  ax.set(ylabel="energy (in wave numbers)")
  ax.legend()

def main():
  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8,7.2), dpi=300)
  fig.tight_layout(pad=5.0)
  
  plotCO(ax1)
  plotN2(ax2)
  
  plt.suptitle("Comparison of potential approximations", fontweight='bold')
  
  plt.savefig("plot.png")
  
if __name__ == "__main__":
  main()
