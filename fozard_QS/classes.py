import numpy as np
# A collection of all the classes currently.

### INITIALIZE (Fozard)
vortex_length = vl = 10 # Length of lattice (micro m), square
Lx, Ly, Lz = 10*vl, 10*vl, 100*vl # Length of Biofilm
dt = 1/600 # min = 0.1 sec. Only one type of time step unlike Fozard
max_mass_particle = 14700
avg_mass_cell = 410
density_cell = 290
density_eps = 290

max_volume_fraction_cell = 0.52
diffusion_substrate = 40680
diffusion_qsm = diffusion_qsi = 33300
steps_substrate = 12 #Ignored
steps_qsm = steps_qsi = 10 #Ignored

max_mass_eps = density_eps/density_cell * max_mass_particle


# Equation
half_saturation = Ks = 1 # Value?
max_substrate_uptake_rate = Vmax = 1 # Value?
substrate_uptake_rate = Vmax * conc_subst_vortex / (Ks + conc_subst_vortex) * mass_cell



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
                    vortex_arr.append(Vortex(x, y, z, vortex_length))
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
        for vortex in vortex_arr:
            vortex.update()


class Vortex:
    # Square compartments containing particles (Cells, EPS)
    def __init__(self, x, y, z, particle_arr=[], eps_amount=0):
        self.x, self.y, self.z = x, y, z #position
        self.particle_arr = particle_arr
        self.eps_amount = eps_amount
        self.conc_subst = 0 #substrate
        self.conc_qsm = 0 #quorom sensing molecule
        self.conc_qsi = 0 #quorom sensing inhibitor

    def get_pos(self):
        return [self.x, self.y, self.z]

    def update(self):
        for particle in particle_arr:
            particle.update()



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

def model_concentration(conc, conc_neigh_arr, production, diffusion_const):
    # return derivative of concentration.
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


#x = Biofilm()
#y = x.get_vortex(2,3,4)
#z = x.get_vortex_neighbours(y)
#print(y.get_pos(),[r.get_pos() for r in z])
