import numpy as np
import matplotlib.pyplot as plt

def main():
  # For C-O bond
  D_e = 40100 # 1/nm
  r_0 = 1.229 # nm
  a = 2.32 # nm
  k_r = 570.0 # kcal /(mol A^2)
  
  r = np.arange(0, 10, 0.1)
  V_harmonic = 0.5 * k_r * (r - r_0)**2
  V_morse = D_e * (1 - np.exp(-a * (r - r_0)))**2
  
  plt.plot(r, V_morse)
  plt.plot(r, V_harmonic)
  
  
  plt.title("Comparison of harmonic potential energy function vs. Morse potential energy function")
  # plt.xlabel()
  
  plt.ylim(0, 50000)
  plt.show()
  plt.savefig("plot.png")
  
if __name__ == "__main__":
  main()
