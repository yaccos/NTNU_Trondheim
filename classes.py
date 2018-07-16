# Standard Packages
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

# Files
from input_file_fozard import *

# A collection of all the classes currently.



### CLASSES
class Biofilm:
    # Collection of vortexes of biofilm "lattice"
    def __init__(self):
        self.length = [Lx, Ly, Lz] # Total size of biofilm lattice
        self.vortex_arr = self.init_vortex() # Collection of vortexes
        self.time_step = 0


    def init_vortex(self):
        vortex_arr = []
        [nx, ny, nz] = [ int(self.length[i] / vortex_length) for i in range(3)]
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    vortex_arr.append(Vortex(x, y, z, conc_subst = conc_bulk))
        return vortex_arr


    def get_vortex(self, x, y, z):
        # Gets vortex from vortex_arr at pos = x, y, z
        [nx, ny, nz] = [ int(self.length[i] / vortex_length) for i in range(3)]
        if 0 <= z < nz:
            # Continuous boundary condition (x & y direction)
            if x == -1:
                x = nx - 1
            elif x == nx:
                x = 0
            if y == -1:
                y = ny - 1
            elif y == ny:
                y = 0
            index = nx * ny * z + nx * y + x
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


    def update(self):
        #Update all vortexes, which then also updates particles
        self.time_step += 1
        for vortex in self.vortex_arr:
            neigh = self.get_vortex_neighbours(vortex)
            vortex.update(neigh)
        self.bulk_vortex()

    def bulk_vortex(self):
        Nz = int(self.length[2] / vortex_length)
        for vortex in self.vortex_arr:
            if vortex.z == Nz - 1:
                vortex.conc_subst = conc_bulk
                vortex.particle_arr = []
                vortex.pressure = 0

class Vortex:
    # Square compartments containing particles (Cells, EPS)
    def __init__(self, x, y, z, particle_arr=[], eps_amount=0, conc_subst=0, conc_qsm=0, conc_qsi=0):
        self.x, self.y, self.z = x, y, z #position
        self.particle_arr = particle_arr.copy()
        self.eps_amount = eps_amount
        self.conc_subst = self.cs1 = conc_subst #substrate
        self.conc_qsm = self.cqsm1 = conc_qsm #quorom sensing molecule
        self.conc_qsi = self.cqsi1 = conc_qsi #quorom sensing inhibitor
        self.N_p = len(particle_arr)
        self.pressure = self.N_p / (max_particles - self.N_p )


    def get_pos(self):
        return [self.x, self.y, self.z]

    def get_mass(self):
        m = 0
        for particle in self.particle_arr:
            m += particle.mass
        return m

    def update(self, neighbours):
        # Time step. Ignore production currently, to be fixed
        # Creates variables cs1, cqsm1, cqsi1 which stores the new concentrations

        self.conc_subst, self.conc_qsm, self.conc_qsi = self.cs1, self.cqsm1, self.cqsi1

        self._update_particles()
        self._update_eps()
        self._update_displacement(neighbours)
        self._update_concentration(neighbours)

    def add_particle(self, particle):
        self.particle_arr.append(particle)
        N_p = len(self.particle_arr)
        self.pressure = N_p / (max_particles - N_p )


    def _update_particles(self):
        #Particle and prod_subst
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell):
                if particle.mass > max_mass_cell:
                    randf = 0.4 + 0.2*np.random.random() # 0.4 - 0.6
                    
                    new_particle = Particle_Cell( (1-randf) * particle.mass) 
                    particle.set_mass(particle.mass * randf)
                    
                    # Distribute up cells between the particles (depending on num of cells)
                    n = particle.num_up
                    for _ in range(n):
                        tot_num_down = particle.num_down + new_particle.num_down
                        if particle.num_down / tot_num_down < np.random.random():
                            particle.create_up()
                        else:
                            new_particle.create_up()

                    self.add_particle(new_particle)

                v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * particle.mass
                particle.update(self.conc_subst, v, self.conc_qsm, self.conc_qsi)
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
                r = np.random.random()
                for i, p in enumerate(probability):
                    if r < p:
                        index = np.random.randint(Np) # Random index to particles
                        particle = self.particle_arr.pop(index) # Choose random particle
                        Np -= 1
                        neighbours[i].add_particle(particle)
                        break


    def _update_concentration(self, neighbours):
        # Neighbours and concentrations
        cs0, cqsm0, cqsi0 = self.conc_subst, self.conc_qsm, self.conc_qsi  
        prod_subst = 0
        prod_qsm = 0
        prod_qsi = 0 # No production from bacteria

        cs_neigh, cqsm_neigh, cqsi_neigh = [], [], []
        
        for vortex in neighbours:
            cs_neigh.append(vortex.conc_subst)
            cqsm_neigh.append(vortex.conc_qsm)
            cqsi_neigh.append(vortex.conc_qsi)
        
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell): 
                v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * particle.mass
                prod_subst -= v

                u = particle.num_up; d = particle.num_down
                if not cqsm0 == 0: 
                    prod_qsm = Zqu * u * cqsm0 / (Kq + cqsm0) + Zqd * d

        t = (0, dt)
        # self.cs1    = cs0    + dt*model_concentration(cs0,   cs_neigh,   diffusion_subst, prod_subst)
        

        
        arr = odeint(model_concentration, cs0, t, args=(cs_neigh,diffusion_subst, prod_subst) ) 
        self.cs1 = arr[1][0]
        arr =odeint(model_concentration, cqsm0, t, args=(cs_neigh,diffusion_qsm, prod_qsm) ) 
        self.cqsm1 = arr[1][0]
        arr =odeint(model_concentration, cqsi0, t, args=(cs_neigh,diffusion_qsi, prod_qsi) ) 
        self.cqsi1 = arr[1][0]
        
 

class Particle_Cell:
    def __init__(self, mass):
        self.mass = mass #changes over time
        self.num_down = math.ceil( mass / avg_mass_cell) #Down-regulated cell
        self.num_up = 0 #Up regulated cell (produces more EPS)
        
    def update(self, conc_subst, v, conc_qsm, conc_qsi):
        # Model eating of substrate and switching states (stochastic)
        self.set_mass(self.mass + dt * model_cell_mass(conc_subst, self.mass, v) )

        for _ in range(self.num_down):
            randf = np.random.random()
            if randf < dt * probability_down2up(conc_qsm, conc_qsi):
                self.create_up()

        for _ in range(self.num_up):
            randf = np.random.random()
            if randf < dt * probability_up2down(conc_qsm, conc_qsi):
                self.create_down()



    def get_cells(self):
        #Returns the number of cells in particle
        return self.num_down + self.num_up

    def set_mass(self, mass):
        # Creates down regulated cell automatically. 
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
    # Not implemented
    def __init__(self, mass):
        self.mass = mass
        return 0

    def update(self):
        # No changes are made after creation
        return 



### MODELS (Time derivatives)

def model_concentration(conc, t, conc_neigh_arr, diffusion_const, production = 0):
    # return derivative of concentration
    # diffusion_const = diffusion_subst, diffusion_qsi, diffusion_qsm
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


def probability_down2up(conc_qsm, conc_qsi):
    a = alpha
    c = gamma
    qsm = conc_qsm
    qsi = conc_qsi

    Qplus = alpha * qsm / (1 + c * (qsm + qsi) )
    return Qplus


def probability_up2down(conc_qsm, conc_qsi):
    b = beta
    c = gamma
    qsm = conc_qsm
    qsi = conc_qsi

    Qminus = b * (1 + gamma*qsi) / (1 + gamma*(qsm + qsi) )
    return Qminus


def time_step(N_times, biofilm):
    # Does N time-steps. Includes a loading "bar".
    loading_bar = [j for j in range(110)]
    for i in range(N_times):
        biofilm.update()
        if i*100.0/N_times >= loading_bar[0]:
            x = loading_bar.pop(0)
            print("%i %%" % x, end='\t\t')
            estimate_time(biofilm, N_times-i)

from time import time
def estimate_time(bf, N):
    # Calculate estimated time (using 1 iterations)
    start_time = time()
    bf.update()
    est_time = N * (time() - start_time)
    print("Estimated time left: %.2f h = %.0f min" % (est_time/3600, est_time/60)  )
    


