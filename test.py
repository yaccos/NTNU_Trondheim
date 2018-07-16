from input_file_fozard import *
from print_file import print_fozard
import classes
from numpy.random import random
from time import time
import datetime
now = datetime.datetime.now()

    
    

### INITIALIZE BIOFILM
bf = classes.Biofilm()
for vortex in bf.vortex_arr:
    if vortex.z == 0: # 1-2 cells in 10 particles per vortex
        vortex.particle_arr = [classes.Particle_Cell(400 + 400*random()) for _ in range(10) ]

print_fozard(bf, "data/before_small.dat", 0)
print("Start-time:", now.strftime("%Y-%m-%d %H:%M"))

### CALCULATE N STEPS
start_time = time()
classes.time_step(N, bf)
delta_time = time() - start_time

# print
print("Model time (min): %.2f" % (bf.time_step * dt) )
print("Used time(min): %.2f" % (delta_time/60) )
print_fozard(bf, "data/small.dat", delta_time)
