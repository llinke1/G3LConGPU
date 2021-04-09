## @package combineGtilde.py
# Python script that reads in Gtilde calculated for different tiles
# and combines them
# Can also give out jackknife samples
#
# @param argv[0] Jackknife sample that is left out, 0 if all tiles are taken
# @param argv[1] Number of tiles
# @param argv[2] Outputfilename for Gtilde
# @param argv[3:N+3] Gtilde files
#
# @author Laila Linke, llinke@astro.uni-bonn.de


import sys
import numpy as np

# Parse cmd line
jn=int(sys.argv[1])
N=int(sys.argv[2])
fn_out=sys.argv[3]
fns_gtilde=sys.argv[4:N+4]
# Take out jackknife
if(jn > 0):
    fns_gtilde=fns_gtilde.pop(jn-1)

# Combine Gtilde
# Read in first
data=np.loadtxt(fns_gtilde[0])
print("Reading ", fns_gtilde[0])
theta1=data[:,0]
theta2=data[:,1]
phi=data[:,2]
dtheta1=data[:,3]
dtheta2=data[:,4]
dphi=data[:,5]
G_real=data[:,6]
G_imag=data[:,7]
weight=data[:,8]

for fn_gtilde in fns_gtilde[1:]:
    print("Reading ", fn_gtilde)
    data=np.loadtxt(fn_gtilde)
    newG_real=data[:,6]
    newG_imag=data[:,7]
    newWeight=data[:,8]
    for i, G in enumerate(G_real):
        w=weight[i]+newWeight[i]
        if(w!=0):
            G_real[i] = G_real[i]*weight[i]/w + newG_real[i]*newWeight[i]/w
            G_imag[i] = G_imag[i]*weight[i]/w + newG_imag[i]*newWeight[i]/w
            
            weight[i] = w

G_real=np.nan_to_num(G_real)
G_imag=np.nan_to_num(G_imag)

# Output
result=np.column_stack((theta1, theta2, phi, dtheta1, dtheta2, dphi, G_real, G_imag, weight))
np.savetxt(fn_out, result)
