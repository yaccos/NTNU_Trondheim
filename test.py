from input_file_fozard import *
from print_file import print_fozard
from classes import Biofilm, Vortex, Particle_Cell, Particle_EPS
from np.random import random

print_fozard()

### Tests
bf = Biofilm()



def fozard_init():
    return



# INITIALIZE BIOFILM
for vortex in bf.vortex_arr:
    vortex.conc_subst = vortex.cs1 = conc_bulk
    vortex.conc_qsm = 0.0
    vortex.conc_qsi = 0.0
    vortex.eps_amount = 0.0
    if vortex.z == 0: # 1-2 cells in 10 particles per vortex
        vortex.particle_arr = [Particle_Cell(400 + 400*random()) for _ in range(10) ]

time_step(N, bf)

print("time (min):", bf.time_step * dt )

