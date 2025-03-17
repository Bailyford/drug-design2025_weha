import numpy as np
import sys
energy = np.loadtxt("energy.dat", usecols=[8], skiprows=0)
np.savetxt(sys.stdout, energy)

