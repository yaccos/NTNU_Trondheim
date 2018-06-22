import numpy as np
# A collection of all the classes currently.

# Initialize
vortex_length = vl = 10 # Length of lattice (micro m), square
Lx, Ly, Lz = 10*vl, 10*vl, 100*vl # Length of Biofilm
dt = 1/600 # min = 0.1 sec. Only one type of time step
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
        #Update all vortexes -> particles
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
        self.num_cells = math.ceil( mass / avg_mass) #changes with mass
        self.num_down = self.num_cells
        self.num_up = 0
        

class Particle_EPS:
    # Not implemented
    def __init__(self):
        return 0


def model_subst(conc, conc_neigh_arr, production):
    # return derivative of conc_subst
    Ds = diffusion_substrate
    l = vortex_length
    f = production
    n = len(conc_neigh_arr)

    dcdt = Ds/l**2 * (sum(conc_neigh_arr) - n*conc) + production/l**3
    return dcdt

def model_qsm(conc, conc_neigh_arr, production):
    # return derivative of conc_subst
    Ds = diffusion_substrate
    l = vortex_length
    f = production
    n = len(conc_neigh_arr)

    dcdt = Ds/l**2 * (sum(conc_neigh_arr) - n*conc) + production/l**3
    return dcdt

def model_qsi(conc, conc_neigh_arr, production):
    # return derivative of conc_subst
    Ds = diffusion_substrate
    l = vortex_length
    f = production
    n = len(conc_neigh_arr)

    dcdt = Ds/l**2 * (sum(conc_neigh_arr) - n*conc) + production/l**3
    return dcdt

#x = Biofilm()
#y = x.get_vortex(2,3,4)
#z = x.get_vortex_neighbours(y)
#print(y.get_pos(),[r.get_pos() for r in z])
