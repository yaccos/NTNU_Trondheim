import numpy as np
# A collection of all the classes currently.

### INITIALIZE (Fozard)
vortex_length = vl = 10 # Length of lattice (micro m), square
Lx, Ly, Lz = 2*vl, 2*vl, 2*vl # Length of Biofilm
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
half_saturation = Ks = 1 # Value?
max_substrate_uptake_rate = Vmax = 1 # Value?
#substrate_uptake_rate = Vmax * conc_subst_vortex / (Ks + conc_subst_vortex) * mass_cell



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
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    vortex_arr.append(Vortex(x, y, z))
        return vortex_arr


    def get_vortex(self, x, y, z):
        # Inefficient.
        # Gets vortex from vortex_arr at pos = x, y, z
        for vortex in self.vortex_arr:
            [vx, vy, vz] = vortex.get_pos()
            if x == vx and y == vy and z == vz:
                return vortex
        #print("Vortex at ", x, y, z, "not found")#For debug
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
        self.time_step += 1
        #Update all vortexes, which then also updates particles
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
        prod_subst, prod_qsi, prod_qsm = 0, 0, 0

        cs0, cqsm0, cqsi0 = self.cs1, self.cqsm1, self.cqsi1
        cs_neigh, cqsm_neigh, cqsi_neigh = [], [], []

        for vortex in neighbours:
            print
            cs_neigh.append(vortex.conc_subst)
            cqsm_neigh.append(vortex.conc_qsm)
            cqsi_neigh.append(vortex.conc_qsi)

        cs1 = cs0 + dt*model_concentration(cs0, cs_neigh, diffusion_subst, prod_subst)
        cqsm1 = cqsm0 + dt*model_concentration(cqsm0, cqsm_neigh, diffusion_qsm, prod_qsi)
        cqsi1 = cqsi0 + dt*model_concentration(cqsi0, cqsi_neigh, diffusion_qsi, prod_qsm)
        for particle in self.particle_arr:
            particle.update()

        self.cs1 = cs1        
        self.cqsm1 = cqsm1
        self.cqsi1 = cqsi1
        

class Particle_Cell:
    def __init__(self, mass):
        self.mass = mass #changes over time
        self.num_down = math.ceil( mass / avg_mass) #Down-regulated cell
        self.num_up = 0 #Up regulated cell (produces more EPS)
        
    def update(self):
        return 0

    def get_cells(self):
        #Returns the number of cells in particle
        return self.num_down + self.num_up


class Particle_EPS:
    # Not implemented
    def __init__(self, mass):
        self.mass = mass
        return 0

    def update(self):
        return 0


### MODELS (Time derivatives)

def model_concentration(conc, conc_neigh_arr, diffusion_const, production = 0):
    # return derivative of concentration. Temporarily production = 0
    # diffusion_const = diffusion_subst, diffusion_qsi, diffusion_qsm
    D = diffusion_const
    l = vortex_length
    f = production
    n = len(conc_neigh_arr)
    dcdt = D/l**2 * (sum(conc_neigh_arr) - n*conc) + production/l**3
    return dcdt


def model_particle_mass(Ymax, substrate_uptake_rate, maintenance_rate, mass_particle):
    # Yield=ymax
    # Returns the derivative of particle mass
    v = substrate_uptake_rate
    m = maintenance_rate
    M = mass_particle

    dMdt = Ymax * (v - m*M)
    return dMdt


bf = Biofilm()

vortex = bf.vortex_arr[0]
vortex.conc_subst = vortex.cs1 = 1000

print("Before")
for vortex in bf.vortex_arr:
    print(vortex.conc_subst)
bf.update()
print("After")
for vortex in bf.vortex_arr:
    print(vortex.cs1)


#y = x.get_vortex(2,3,4)
#z = x.get_vortex_neighbours(y)
#print(y.get_pos(),[r.get_pos() for r in z])
