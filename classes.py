# Standard Packages
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

# Files
from input_file import * #If there's a variable you cannot find, it's probably here.

### A collection of all the classes currently:
# Biofilm, collection of vortexes
# Vortex, "Boxes" containing particles and concentrations
# Particle_Cell, Particles of "n" cells, up/down regulated
# Particle_EPS, Biofilm particles



### CLASSES
class Biofilm:
    # Collection of vortexes of biofilm "lattice"
    # TODO Concentration array instead of Vortexes
    def __init__(self):
        self.length = [Lx, Ly, Lz] # Total size of biofilm lattice
        self.vortex_arr = self._init_vortex() # Collection of vortexes in 1D list
        self.time_step = 0


    def update(self):
        # Update all vortexes, which then also updates particles
        self.time_step += 1
        [Nx, Ny, Nz] = self.get_box_size()
        for vortex in self.vortex_arr:
            neigh = self.get_vortex_neighbours(vortex)
            vortex.update(neigh, Nz)


    def get_box_size(self):
        # Returns the number of boxes in each direction [Nx, Ny, Nz]
        return [int(self.length[i] / vortex_length) for i in range(3) ]


    def get_vortex(self, x, y, z):
        # Gets vortex from vortex_arr at pos = x, y, z
        [Nx, Ny, Nz] = self.get_box_size()
        if 0 <= z < Nz:
            x, y = self._cbc(x,y)
            index = Nx * Ny * z + Nx * y + x
            return self.vortex_arr[index]
        else:
            return None

    def get_vortex_neighbours(self, vortex):
        # Neighbours when sides are touching, not when only corners
        # Returns an array of all the neighbour vortexes of the vortex input.
        neighbours = []
        [x,y,z] = vortex.get_pos()
        for pos in [[x+1,y,z], [x-1,y,z], [x,y+1,z], [x,y-1,z], [x,y,z+1], [x,y,z-1]]:
            xv, yv, zv = pos[0], pos[1], pos[2]
            vortex = self.get_vortex(xv,yv,zv)
            if vortex != None:
                neighbours.append( vortex )
        return neighbours


    def _init_vortex(self):
        # Creates vortex at every xyz with substrate concentration of the bulk liquid
        [Nx, Ny, Nz] = self.get_box_size()
        vortex_arr = []
        for z in range(Nz):
            for y in range(Ny):
                for x in range(Nx):
                    vortex_arr.append(Vortex(x, y, z, conc_subst=conc_bulk))
        return vortex_arr


    def _cbc(self,x,y):
        # Continuous boundary condition (x & y direction)
        [Nx, Ny, Nz] = self.get_box_size()
        if x == -1:
            x = Nx - 1
        elif x == Nx:
            x = 0
        if y == -1:
            y = Ny - 1
        elif y == Ny:
            y = 0
        return x,y
            


class Vortex:
    # Square compartments containing particles (Cells, EPS)
    def __init__(self,
            x, y, z, 
            cell_mass_arr=[], eps_mass_arr=[],
            eps_amount=0, conc_subst=0, conc_qsm=0):

        self.x, self.y, self.z = x, y, z #position
        
        # Same index refers to same particle
        self.cell_mass_nparr = np.array(cell_mass_arr)
        self.cell_up_arr = [0] * len(cell_mass_arr) 
        self.cell_down_arr =  np.ceil(self.cell_mass_nparr/avg_mass_cell).astype(int).tolist()


        self.eps_mass_arr = eps_mass_arr.copy()
        self.eps_amount = eps_amount
        self.conc_subst = self.cs1 = conc_subst #substrate
        self.conc_qsm = self.cqsm1 = conc_qsm #quorom sensing molecule

        self.pressure = 0
        self.update_pressure()

        self.counter = 0

    def get_num_particles(self):
        return len(self.eps_mass_arr) + len(self.cell_mass_nparr)

    def get_pos(self):
        return [self.x, self.y, self.z]

    def get_mass(self):
        return sum(self.eps_mass_arr) + sum(self.cell_mass_nparr)

    def update(self, neighbours, Nz):
        self.counter += 1
        # Time step.
        # Creates variables cs1, cqsm1 which stores the new concentrations
        # TODO Is the concentration before/after correctly handled?
        if self.z == Nz - 1: # If in bulk-liquid
            self.conc_subst = conc_bulk
            self.conc_qsm = 0
            self.pressure = 0

        self.conc_subst, self.conc_qsm = self.cs1, self.cqsm1 

        self._update_particles()
        self._update_eps()
        self._update_displacement(neighbours)
        self._update_concentration(neighbours)

    def add_eps(self, mass):
        self.eps_mass_arr.append(mass)
        self.update_pressure()

    
    def add_cell(self, cell_mass, down, up):
        self.cell_mass_nparr = np.append(self.cell_mass_nparr, cell_mass)
        self.cell_up_arr.append(up)
        self.cell_down_arr.append(down)
        self.update_pressure()


    def update_pressure(self):
        N = self.get_num_particles()
        self.pressure = N / (max_particles - N )

    def create_up(self,index):
        i=index
        self.cell_up_arr[i] += 1
        self.cell_down_arr[i] -= 1
    def create_down(self,index):
        i=index
        self.cell_up_arr[i] -= 1
        self.cell_down_arr[i] += 1

    def _update_particles(self):
        #Particle and prod_subst
        # Currently the slowest part
        #TODO Use np.random.binomial(n, probability) instead of for loop
        #TODO Use eps and cell array with numpy/scipy
        index_split_mass = [i for i in range(len(self.cell_mass_nparr)) if self.cell_mass_nparr[i] > max_mass_cell]

        for i in index_split_mass:
            randf = 0.4 + 0.2* random() # 0.4 - 0.6
            mass = self.cell_mass_nparr[i]
            self.add_cell(mass * (1-ranf), ceil(mass/avg_mass_cell), 0)
            self.cell_mass_nparr[i] *= randf

            up, down = self.cell_up_arr[i], self.cell_down_arr[i]
            transfer_up2new_cell = np.random.hypergeometric(up, down, up)
            for _ in range(transfer_up2new_cell):
                self.create_down(i)
                self.create_up(-1) # -1 = new cell

        self._update_cell()
         

    def _update_cell(self):
        # TODO This function is from particle_cell
        # Model eating of substrate and switching states (stochastic)
        self._update_mass()
        pd2u = probability_down2up(self.conc_qsm)*dt
        pu2d = probability_up2down(self.conc_qsm)*dt 
        sucess2up = np.random.binomial(self.cell_up_arr, pd2u)
        sucess2down = np.random.binomial(self.cell_down_arr, pu2d)

        for i in range(len(sucess2up)):  
            for _ in range(int(sucess2up[i])):
                self.create_up(i)
        for i in range(len(sucess2down)):  
            for _ in range(int(sucess2down[i])):
                self.create_down(i)



    def _update_mass(self):
        v = Vmax * self.conc_subst / (Ks + self.conc_subst) * self.cell_mass_nparr #Substrate uptake rate, np array
        # Updates mass according to model
        new_mass_nparr = self.cell_mass_nparr + dt * model_cell_mass(self.conc_subst, self.cell_mass_nparr, v) 
        self.set_mass(new_mass_nparr)
        

    def set_mass(self,mass_nparr):
        # Updates cell count when mass updates
        self.cell_mass_nparr = mass_nparr
        upanddown_mass_arr = avg_mass_cell *(np.array(self.cell_up_arr) + np.array(self.cell_down_arr))
        
        x = (mass_nparr - avg_mass_cell > upanddown_mass_arr) #Points where mass is high enough to warrant another cell
        for i in x.astype(int):
            self.cell_down_arr[i] += 1
        
        #Prioritize removing down regulated
        y =(mass_nparr+avg_mass_cell < upanddown_mass_arr )
        for i in x.astype(int):
            if self.cell_down_arr[i] > 0:
                self.cell_down_arr[i] -= 1
            elif self.cell_up_arr[i] > 0:
                self.cell_up_arr[i] -= 1


    def _update_eps(self):
        #EPS production from bacteria
        self.eps_amount += dt * model_eps(self.cell_down_arr, self.cell_up_arr)

        # If mass surpasses an amount, create particle from that mass
        if self.eps_amount > max_mass_eps:
            self.eps_mass_arr.append( max_mass_eps )
            self.eps_amount -= max_mass_eps


    def _update_displacement(self, neighbours):
        # Displacement of particles (Pressure)
        self.update_pressure()
        mu = transfer_coefficient

        delta_Np = 0; tot_diff_pressure = 0
        Np = self.get_num_particles()
        for vortex in neighbours:
            Npp = vortex.get_num_particles()
            if self.pressure > vortex.pressure and Np > Npp:
                delta_Np += math.floor( mu* (self.pressure - vortex.pressure) * (Np - Npp))
            tot_diff_pressure += self.pressure - vortex.pressure

        probability = np.zeros(len(neighbours))
        for i, vortex in enumerate(neighbours):
            if self.pressure <= vortex.pressure:
                probability[i] = 0
            else:
                probability[i] = (self.pressure - vortex.pressure) / tot_diff_pressure

        # From discrete distribution to cumulative distribution
        for i in range(len(probability) - 1):
            probability[i+1] += probability[i] # [0.1, 0.2, 0.4, 0.3] -> [0.1, 0.3, 0.7, 1] 
       

        for _ in range(delta_Np):
            # Choose either cell or EPS particle
            if np.random.hypergeometric(len(self.cell_mass_nparr), len(self.eps_mass_arr), 1): #True -> cell
                for i, p in enumerate(probability):
                    if r <= p:
                        index = np.random.randint(len(self.cell_mass_nparr)) # Random index to particles
                        mass = self.cell_mass_nparr[index] 
                        self.cell_mass_nparr = np.delete(self.cell_mass_nparr, index)
                        up = self.cell_up_arr.pop[index] 
                        down = self.cell_down_arr.pop[index]

                        neighbours[i].add_cell(mass, down, up)
                        break
            else:
                for i, p in enumerate(probability):
                    if r <= p:
                        x = self.eps_mass_arr.pop()
                        neighbours[i].add_eps(x)
                        break


    def _update_concentration(self, neighbours):
        # Initialize
        cs0, cqsm0 = self.conc_subst, self.conc_qsm
        prod_subst = 0
        prod_qsm = 0
        
        # Neighbours 
        cs_neigh, cqsm_neigh = [], []
        for vortex in neighbours:
            cs_neigh.append(vortex.conc_subst)
            cqsm_neigh.append(vortex.conc_qsm)
        
        # Production
        v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * self.cell_mass_nparr
        prod_subst -= np.sum(v)
        u = sum(self.cell_up_arr)
        d = sum(self.cell_down_arr)
        prod_qsm = Zqu * u * cqsm0 / (Kq + cqsm0) + Zqd * d

        # Time step
        t = (0, dt)
        
        arr = odeint(model_concentration, cs0, t, args=(cs_neigh,diffusion_subst, prod_subst) ) 
        self.cs1 = arr[1][0]
        arr =odeint(model_concentration, cqsm0, t, args=(cqsm_neigh,diffusion_qsm, prod_qsm) ) 
        self.cqsm1 = arr[1][0]
        

### MODELS (Time derivatives)

def model_concentration(conc, t, conc_neigh_arr, diffusion_const, production=0):
    # return derivative of concentration
    # diffusion_const = diffusion_subst, diffusion_qsm
    D = diffusion_const
    l = vortex_length
    f = production
    n = len(conc_neigh_arr)

    dcdt = D/l**2 * (sum(conc_neigh_arr) - n*conc) + production/l**3
    return dcdt


def model_cell_mass(conc_subst, mass, v):
    # Returns the derivative of particle mass
    m = maintenance_rate
    Y = max_yield
    M = mass

    dMdt = Y * (v - m*M)
    return dMdt


def model_eps(cell_up_arr, cell_down_arr):
    # Production of eps from cells (different prod for up/down given by Zd, Zu)
    dEdt = Zed * sum(cell_down_arr) + Zeu * sum(cell_up_arr)

    return dEdt


def probability_down2up(conc_qsm):
    a = alpha
    y = gamma
    qsm = conc_qsm

    Qplus = alpha * qsm / (1 + y * qsm)
    return Qplus


def probability_up2down(conc_qsm):
    b = beta
    y = gamma
    qsm = conc_qsm

    Qminus = b / (1 + y*qsm )
    return Qminus

### OTHER

import print_file
def time_step(N_times, biofilm):
    # Does N time-steps. Includes a loading "bar".
    loading_bar = [j for j in range(110)]
    for i in range(N_times):
        biofilm.update()
        if i*100.0/N_times >= loading_bar[0]:
            x = loading_bar.pop(0)
            print("%i %%" % x, end='\t\t')
            estimate_time(biofilm, N_times-i)
            # Every 10% print to a file
            if int(x) % 10 == 0:
                print_file.print_output("data/temporary" + str(x) + ".dat", biofilm)


from time import time
def estimate_time(bf, N):
    # Calculate estimated time (using 1 iterations)
    start_time = time()
    bf.update()
    est_time = N * (time() - start_time)
    hour = int(est_time/3600)
    minute = int(est_time/60 - hour*60) 
    print("ETA: %i h  %i min" % (hour, minute) )

