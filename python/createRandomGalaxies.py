import numpy as np
import sys

tile_length=float(sys.argv[1]) #arcmin
outfn=sys.argv[2]

num_density=30 # arcmin^-2

N=int(tile_length*tile_length*num_density) #Number of galaxies


x=np.random.uniform(size=N)*tile_length
y=np.random.uniform(size=N)*tile_length
z=np.random.uniform(size=N)*2

result=np.column_stack((x, y, z))
np.savetxt(outfn, result)
