import numpy as np
import matplotlib.pyplot as plt
import math
# A collection of all the classes currently.


### INITIALIZE (Fozard)
vortex_length = vl = 10 # Length of lattice (micro m), square
Lx, Ly, Lz = 20*vl, 10*vl, 1*vl # Length of Biofilm
dt = 1/60000 # min = 0.001 sec. Only one type of time step unlike Fozard
max_mass_particle = 14700
avg_mass_cell = 410
density_cell = 290
density_eps = 290

max_volume_fraction_cell = 0.52
diffusion_subst = 40680
diffusion_qsm = diffusion_qsi = 33300
steps_substrate = 12 #Ignored
steps_qsm = steps_qsi = 10 #Ignored

max_mass_eps = density_eps/density_cell * max_mass_particle


# Equation
half_saturation = Ks = 2.34e-3 # Value?
max_substrate_uptake_rate = Vmax = 0.046 # Value?
max_yield = 0.444 # Value?
maintenance_rate = 6e-4

Zd = 1e-6   # Down-regulated EPS production
Zu = 1e-3   # Up-regulated EPS production



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
        self.particle_arr = particle_arr
        self.eps_amount = eps_amount
        self.conc_subst = self.cs1 = conc_subst #substrate
        self.conc_qsm = self.cqsm1 = conc_qsm #quorom sensing molecule
        self.conc_qsi = self.cqsi1 = conc_qsi #quorom sensing inhibitor


    def get_pos(self):
        return [self.x, self.y, self.z]


    def update(self, neighbours):
        # Time step. Ignore production currently, to be fixed
        # Creates variables cs1, cqsm1, cqsi1 which stores the new concentrations
        prod_subst = 0 # This should work
        prod_qsi, prod_qsm = 0, 0, 0 #Not implemented yet

        self.conc_subst, self.conc_qsm, self.conc_qsi = cs0, cqsm0, cqsi0 = self.cs1, self.cqsm1, self.cqsi1

        #Particle and prod_subst
        for particle in self.particle_arr:
            if isinstance(particle, Particle_Cell):
                self.eps_amount += dt * (Zd * particle.num_down + Zu * particle.num_up)
                if particle.mass > max_mass_cell:
                    randf = 0.4 + 0.2*np.random.random() # 0.4 - 0.6
                    
                    new_particle = Particle_Cell((1-ranf) * particle.mass)) 
                    particle.set_mass(particle.mass * randf)

                    for i in range(

                    particle_arr.append(new_particle)

            v = substrate_uptake_rate = Vmax * conc_subst / (Ks + conc_subst) * mass
            particle.update(self.conc_subst, v)
            prod_subst -= v


        # If mass surpasses an amount, create particle from that mass
        if self.eps_amount > max_mass_eps:
            self.particle_arr.append( Particle_EPS(max_mass_eps) )
            self.eps_amount -= max_mass_eps


        # Neighbours and concentrations
        cs_neigh, cqsm_neigh, cqsi_neigh = [], [], []
        for vortex in neighbours:
            cs_neigh.append(vortex.conc_subst)
            cqsm_neigh.append(vortex.conc_qsm)
            cqsi_neigh.append(vortex.conc_qsi)

        self.cs1 = cs0 + dt*model_concentration(cs0, cs_neigh, diffusion_subst, prod_subst)
        self.cqsm1 = cqsm0 + dt*model_concentration(cqsm0, cqsm_neigh, diffusion_qsm, prod_qsi)
        self.cqsi1 = cqsi0 + dt*model_concentration(cqsi0, cqsi_neigh, diffusion_qsi, prod_qsm)



class Particle_Cell:
    def __init__(self, mass):
        self.mass = mass #changes over time
        self.num_down = math.ceil( mass / avg_mass_cell) #Down-regulated cell
        self.num_up = 0 #Up regulated cell (produces more EPS)
        
    def update(self, conc_subst):
        # Model eating of substrate
        self.set_mass(self.mass + dt * model_cell_mass(conc_subst, self.mass) )

    def get_cells(self):
        #Returns the number of cells in particle
        return self.num_down + self.num_up

    def set_mass(self, mass):
        # Creates down regulated cell automatically. Should be stochastic
        self.mass = mass
        while mass-avg_mass_cell > avg_mass_cell * self.get_cells(): 
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
                z_plane[vortex.x, vortex.y] += particle.mass
    plt.imshow(z_plane)
    plt.show()
 
### Tests
bf = Biofilm()

vortex = bf.vortex_arr[44]
vortex.conc_subst = vortex.cs1 = 1
vortex.particle_arr = [Particle_Cell(1)]

plot2d_cellmass(bf)
for i in range(100):
    bf.update()
plot2d_cellmass(bf)
print(vortex.particle_arr[0].mass)

