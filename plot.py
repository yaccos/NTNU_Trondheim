import matplotlib.pyplot as plt
from input_file import vortex_length
import numpy as np

def plot2d_subst(ax, biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Nx, Ny, Nz = biofilm.grid_size
    z_plane = np.zeros(Nx, Ny)
    for x in range(Nx):
        for y in range(Ny):
            z_plane[x, y] = biofilm.conc_subst[x, y, 0]
    return ax.imshow(z_plane)


def plot2d_cellmass(ax, biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Nx, Ny, Nz = biofilm.length
    z_plane = np.zeros([Nx, Ny])
    for x in range(Nx):
        for y in range(Ny):
            z_plane[x, y] += biofilm.total_cell_mass[x, y, 0]
    return ax.imshow(z_plane)

 
def plot3d(ax, biofilm):
    min_val = 1e8
    max_val = 1e-8
    Nx, Ny, Nz = biofilm.grid_size
    array = [[[None for z in range(Nz)] for y in range(Ny)] for x in range(Nz)]
    for x in range(Nx):
        for y in range(Ny):
            for z in range(Nz):
                conc = biofilm.conc_subst[[x, y, z]]
                array[x][y][z] = conc
                if max_val < conc:
                    max_val = conc
                if min_val > conc:
                    min_val = conc
    print("Minimum substrate concentration", min_val)
    print("Maximum substrate concentration", max_val)
    for x in range(Nx):
        for y in range(Ny):
            for z in range(Nz):
                if array[x][y][z] < 0:
                    c = '0.0'
                c = '%.1f' % ( 1 - (array[x][y][z]-min_val) / (max_val-min_val))            
                ax.scatter(x, y, z, color=c)
