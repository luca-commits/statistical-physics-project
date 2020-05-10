import numpy as np
import matplotlib.pyplot as plt

def main():
  # For C-O bond, values taken from Cornell (1995) paper (p. 5194)and Morse (1929) paper (p. 62)
  D_e = 91600 # wave no = 1 / cm
  r_0 = 1.229 # Angström
  a = 2.29 # 1 / Angström
  k_r = 12.53e-3 # 1 / (Angström^3)
  
  r = np.arange(0, 5, 0.01)
  
  V_harmonic = 0.5 * k_r * (r - r_0)**2 * 10e8 # wave no = 1 / cm
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2 # wave no = 1 / cm
  
  plt.plot(r, V_morse, label='Morse potential')
  plt.plot(r, V_harmonic, label='Harmonic potential')
  
  plt.ylim(0, 200000)
  
  plt.title("Harmonic potential energy function \n vs. Morse potential energy function for C-O bond")
  plt.xlabel("Nuclear seperation (in Angström)")
  plt.ylabel("Energy (in wave numbers)")
  plt.legend()
  
  plt.gcf().subplots_adjust(left=0.15) # prevent y-Label from being cut off
  
  plt.savefig("plot.png")
  
if __name__ == "__main__":
  main()
