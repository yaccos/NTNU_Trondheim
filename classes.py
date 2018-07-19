# Standard Packages
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

# Files
from input_file import *

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
        for vortex in self.vortex_arr:
            neigh = self.get_vortex_neighbours(vortex)
            vortex.update(neigh)
        self.bulk_vortex()


    def bulk_vortex(self):
        [Nx, Ny, Nz] = self.get_box_size()
        for vortex in self.vortex_arr:
            if vortex.z == Nz - 1:
                vortex.conc_subst = conc_bulk
                vortex.particle_arr = []
                vortex.pressure = 0


    def get_box_size(self):
        # Returns the number of boxes in each direction [Nx, Ny, Nz]
        return [int(self.length[i] / vortex_length) for i in range(3) ]


    def get_vortex(self, x, y, z):
        # Gets vortex from vortex_arr at pos = x, y, z
        [Nx, Ny, Nz] = self.get_box_size()
        if 0 <= z < Nz:
            x, y = _continous_boundary_condition(x,y)
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


    def _continous_boundary_condition(x,y):
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
    def __init__(self, x, y, z, particle_arr=[], eps_amount=0, conc_subst=0, conc_qsm=0):
        self.x, self.y, self.z = x, y, z #position

        #TODO: From particle_arr to eps & cell array
        self.particle_arr = particle_arr.copy()
        

        self.eps_mass_arr = []
        
        # Same index refers to same particle
        self.cell_mass_arr = []
        self.cell_up_arr = []
        self.cell_down_arr = []

        self.eps_amount = eps_amount
        self.conc_subst = self.cs1 = conc_subst #substrate
        self.conc_qsm = self.cqsm1 = conc_qsm #quorom sensing molecule
        self.N_p = len(particle_arr)
        self.pressure = self.N_p / (max_particles - self.N_p )


    def get_pos(self):
        return [self.x, self.y, self.z]

    def get_mass(self):
        m = 0
        for particle in self.particle_arr:
            m += particle.mass
        return m

    def update(self, neighbours, Nz):
        # Time step.
        # Creates variables cs1, cqsm1 which stores the new concentrations
        # Is the concentration before/after correctly handled?
        if self.z == Nz - 1:
            self.conc_subst = 0
            self.conc_qsm = 0
            self.particle_arr = []
            self.pressure = 0


        self.conc_subst, self.conc_qsm = self.cs1, self.cqsm1 

        self._update_particles()
        self._update_eps()
        self._update_displacement(neighbours)
        self._update_concentration(neighbours)

    def add_particle(self, particle):
        self.particle_arr.append(particle)
        self.N_p = len(self.particle_arr)
        self.pressure = self.N_p / (max_particles - self.N_p )


    def _update_particles(self):
        #Particle and prod_subst
        # Currently the slowest part
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell):
                if particle.mass > max_mass_cell:
                    randf = 0.4 + 0.2* random() # 0.4 - 0.6
                    
                    new_particle = Particle_Cell( (1-randf) * particle.mass) 
                    particle.set_mass(particle.mass * randf)
                    
                    # Distribute up cells between the particles (depending on num of cells)
                    u = particle.num_up; d = particle.num_down
                    for _ in range(u):
                        if d == 0:
                            break

                        tot_num_down = d + new_particle.num_down
                        if d / tot_num_down < random():
                            particle.create_up()
                        else:
                            new_particle.create_up()

                    self.add_particle(new_particle)

                v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * particle.mass
                particle.update(self.conc_subst, v, self.conc_qsm)
            else: #isinstance(particle, EPS)
                particle.update() 


    def _update_eps(self):
        #EPS production from bacteria
        cell_arr = []
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell):
                cell_arr.append(particle)
        self.eps_amount += dt * model_eps(cell_arr)

        # If mass surpasses an amount, create particle from that mass
        if self.eps_amount > max_mass_eps:
            self.particle_arr.append( Particle_EPS(max_mass_eps) )
            self.eps_amount -= max_mass_eps


    def _update_displacement(self, neighbours):
        # Displacement of particles (Pressure)
        Np = len(self.particle_arr)
        self.pressure = Np / (max_particles - Np)
        mu = transfer_coefficient

        delta_Np = 0; tot_diff_pressure = 0
        for vortex in neighbours:
            Npp = len(vortex.particle_arr)
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
            if self.particle_arr: #if not empty
                r = random()
                for i, p in enumerate(probability):
                    if r <= p:
                        index = np.random.randint(Np) # Random index to particles
                        particle = self.particle_arr.pop(index) # Choose random particle
                        Np -= 1
                        neighbours[i].add_particle(particle)
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
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell): 
                v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * particle.mass
                prod_subst -= v

                u = particle.num_up; d = particle.num_down
                prod_qsm = Zqu * u * cqsm0 / (Kq + cqsm0) + Zqd * d

        # Time step
        t = (0, dt)
        
        arr = odeint(model_concentration, cs0, t, args=(cs_neigh,diffusion_subst, prod_subst) ) 
        self.cs1 = arr[1][0]
        arr =odeint(model_concentration, cqsm0, t, args=(cqsm_neigh,diffusion_qsm, prod_qsm) ) 
        self.cqsm1 = arr[1][0]
        
 

class Particle_Cell:
    def __init__(self, mass):
        self.mass = mass #changes over time, use set_mass(..)!
        self.num_down = math.ceil( mass / avg_mass_cell) #Down-regulated cell
        self.num_up = 0 #Up regulated cell (produces more EPS)
        
    def update(self, conc_subst, v, conc_qsm):
        # Model eating of substrate and switching states (stochastic)
        _update_mass(conc_subst, v)
        pd2u = probability_down2up(conc_qsm)
        pu2d = probability_up2down(conc_qsm) 
        
        for _ in range(self.num_down):
            if random() < dt * pd2u:
                self.create_up()

        for _ in range(self.num_up):
            if random() < dt * pu2d:
                self.create_down()


    def _update_mass(self, conc_subst, v):
        # Updates mass according to model
        new_mass = self.mass + dt * model_cell_mass(conc_subst, self.mass, v) 
        self.set_mass(new_mass)


    def get_cells(self):
        #Returns the number of cells in particle
        return self.num_down + self.num_up

    def set_mass(self, mass):
        # Creates/Deletes down regulated cell automatically. 
        self.mass = mass
        while mass - avg_mass_cell > avg_mass_cell * self.get_cells(): 
            self.num_down += 1
        while mass+avg_mass_cell < avg_mass_cell * self.get_cells():
            if self.num_down > 0:
                self.num_down -= 1
            elif self.num_up > 0:
                self.num_up -= 1

    def create_up(self):
        self.num_down -= 1
        self.num_up += 1

    def create_down(self):
        self.num_down += 1
        self.num_up -= 1


class Particle_EPS:
    def __init__(self, mass):
        self.mass = mass
        return 0

    def update(self):
        # No changes are made after creation
        return 



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


def model_eps(cell_arr):
    # Production of eps from cells (different prod for up/down given by Zd, Zu)
    dEdt = 0
    for cell in cell_arr:
        d = cell.num_down
        u = cell.num_up

        dEdt += Zed * d + Zeu * u
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

