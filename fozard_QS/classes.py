import numpy as np
import matplotlib.pyplot as plt
import math
# A collection of all the classes currently.


N = int(30e3) # Number of steps
### INITIALIZE (Fozard)
vortex_length = vl = 17 # Length of lattice (micro m), square
Lx, Ly, Lz = 3*vl, 3*vl, 1*vl # Length of Biofilm
dt = 1/6000 # min = 0.01 sec. Only one type of time step unlike Fozard
max_mass_cell = 14700
avg_mass_cell = 410
density_cell = 290
density_eps = 290

max_volume_fraction = 0.52
diffusion_subst = 40680
diffusion_qsm = diffusion_qsi = 33300
steps_substrate = 12 #Ignored
steps_qsm = steps_qsi = 10 #Ignored

max_mass_eps = density_eps/density_cell * max_mass_cell

half_saturation = Ks = 2.34e-3
max_substrate_uptake_rate = Vmax = 0.046
max_yield = 0.444
maintenance_rate = 6e-4

Zd = 1e-6   # Down-regulated EPS production
Zu = 1e-3   # Up-regulated EPS production
max_particles = density_cell * max_volume_fraction * vl**3 / max_mass_cell
print("max particle", math.floor(max_particles))
transfer_coefficient = 0.1



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
                    vortex_arr.append(Vortex(x, y, z))
        return vortex_arr


    def get_vortex(self, x, y, z):
        # Gets vortex from vortex_arr at pos = x, y, z
        [nx, ny, nz] = [ int(self.length[i] / vortex_length) for i in range(3)]
        if 0 <= x < nx and 0 <= y < ny and 0 <= z < nz:       
            index = nx * ny * z + nx * y + x
            return self.vortex_arr[index]
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



class Vortex:
    # Square compartments containing particles (Cells, EPS)
    def __init__(self, x, y, z, particle_arr=[], eps_amount=0, conc_subst=0, conc_qsm=0, conc_qsi=0):
        self.x, self.y, self.z = x, y, z #position
        self.particle_arr = particle_arr.copy()
        self.eps_amount = eps_amount
        self.conc_subst = self.cs1 = conc_subst #substrate
        self.conc_qsm = self.cqsm1 = conc_qsm #quorom sensing molecule
        self.conc_qsi = self.cqsi1 = conc_qsi #quorom sensing inhibitor
        N_p = len(particle_arr)
        self.pressure = N_p / (max_particles - N_p )


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
                self.eps_amount += dt * (Zd * particle.num_down + Zu * particle.num_up)
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

            particle.update(self.conc_subst, v)



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
        prod_qsi, prod_qsm = 0, 0 #Not implemented yet

        cs_neigh, cqsm_neigh, cqsi_neigh = [], [], []
        
        for vortex in neighbours:
            cs_neigh.append(vortex.conc_subst)
            cqsm_neigh.append(vortex.conc_qsm)
            cqsi_neigh.append(vortex.conc_qsi)
        
        for particle in self.particle_arr:
            v = substrate_uptake_rate = Vmax * self.conc_subst / (Ks + self.conc_subst) * particle.mass
            prod_subst -= v

        self.cs1    = cs0    + dt*model_concentration(cs0,   cs_neigh,   diffusion_subst, prod_subst)
        self.cqsm1  = cqsm0  + dt*model_concentration(cqsm0, cqsm_neigh, diffusion_qsm,   prod_qsi)
        self.cqsi1  = cqsi0  + dt*model_concentration(cqsi0, cqsi_neigh, diffusion_qsi,   prod_qsm)
 

class Particle_Cell:
    def __init__(self, mass):
        self.mass = mass #changes over time
        self.num_down = math.ceil( mass / avg_mass_cell) #Down-regulated cell
        self.num_up = 0 #Up regulated cell (produces more EPS)
        
    def update(self, conc_subst, v):
        # Model eating of substrate
        self.set_mass(self.mass + dt * model_cell_mass(conc_subst, self.mass, v) )

    def get_cells(self):
        #Returns the number of cells in particle
        return self.num_down + self.num_up

    def set_mass(self, mass):
        # Creates down regulated cell automatically. Should be stochastic
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



class Particle_EPS:
    # Not implemented
    def __init__(self, mass):
        self.mass = mass
        return 0

    def update(self):
        # No changes are made after creation
        return 



### MODELS (Time derivatives)

def model_concentration(conc, conc_neigh_arr, diffusion_const, production = 0):
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


### PLOT
def plot2d_subst(biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Lx, Ly, Lz = biofilm.length
    z_plane = np.zeros([int(Lx/vortex_length), int(Ly/vortex_length)])
    for vortex in biofilm.vortex_arr:
        if vortex.z == 0:
            z_plane[vortex.x, vortex.y] = vortex.cs1
    plt.imshow(z_plane)
    plt.show()


def plot2d_cellmass(biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Lx, Ly, Lz = biofilm.length
    Nx, Ny = int(Lx/vortex_length), int(Ly/vortex_length)
    z_plane = np.zeros([Nx, Ny])
    for vortex in biofilm.vortex_arr:
        if vortex.z == 0:
            for particle in vortex.particle_arr:
                z_plane[vortex.x, vortex.y] += 1
    plt.imshow(z_plane)
    plt.colorbar()
    plt.show()
 
### Tests
bf = Biofilm()

vortex = bf.vortex_arr[4]
vortex.conc_subst = vortex.cs1 = 100
vortex.particle_arr = [Particle_Cell(1000*avg_mass_cell) for _ in range(20)]
n = bf.get_vortex_neighbours(vortex)
n[0].add_particle(Particle_Cell(100*avg_mass_cell))

bf.update()
print(vortex.particle_arr[0].get_cells())
loading_bar = [N * 0.05 * i for i in range(21)]
for i in range(N):
    bf.update()
    if i == loading_bar[0]:
        loading_bar.pop(0)
        print(int(i/N * 100), "%")
    
print(vortex.particle_arr[0].get_cells())
print("time (min):", bf.time_step * dt )
plot2d_cellmass(bf)

