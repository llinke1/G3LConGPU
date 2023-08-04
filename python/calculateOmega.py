## @package calculateOmega.py
# Python script that reads in paircounts calculated with C code
# and outputs two-point angular correlation function of Galaxies
#
# Paircounts are read in, normalized, and then calculates Omega
# Uses Landy-Szalay for Auto-Correlations and Szalay-Szapudi for
# Cross Correlations
#
# @param argv[1] Boolean, which is 1 for auto- and 0 for cross-corrs
# @param argv[2] N, Number of tiles
# @param argv[3] Outputfilename for omega
# @param argv[4:N+4] Filenames of <DD> paircount
# @param argv[N+4:2N+4] Filenames of <RR> paircount
# @param argv[2N+4:3N+4] Filenames of <D1R2> paircount
# @param argv[3N+4:4N+4] Filenames of <D2R1> paircount (only needed for cross-corrs)
# @author Laila Linke, llinke@astro.uni-bonn.de


import sys
import numpy as np

# Parse cmd line
is_auto=bool(sys.argv[1])
N=int(sys.argv[2])
fn_out=sys.argv[3]
fns_DD=sys.argv[4:N+4]
fns_RR=sys.argv[N+4:2*N+4]
fns_D1R2=sys.argv[2*N+4:3*N+4]

print(fns_DD)

if (not is_auto):
    fns_D2R1=sys.argv[3*N+1:4*N]

# Combine DD

data_DD=np.loadtxt(fns_DD[0])
thetas=data_DD[:,0] # Bin centers for omega
dd=data_DD[:,1] # DD-Paircount
norm=data_DD[0,2] # Normalization

for i, fn in enumerate(fns_DD[1:]):
    data_DD=np.loadtxt(fn) # Read in data_DD
    dd+=data_DD[:,1] # Add paircount
    norm+=data_DD[0,2] # Add normalization

dd/=norm

# Combine RR

data_RR=np.loadtxt(fns_RR[0])
thetas=data_RR[:,0] # Bin centers for omega
rr=data_RR[:,1] # DD-Paircount
norm=data_RR[0,2] # Normalization

for i, fn in enumerate(fns_RR[1:]):
    data_RR=np.loadtxt(fn) # Read in data_RR
    rr+=data_RR[:,1] # Add paircount
    norm+=data_RR[0,2] # Add normalization

rr/=norm

# Combine D1R2

data_D1R2=np.loadtxt(fns_D1R2[0])
thetas=data_D1R2[:,0] # Bin centers for omega
d1r2=data_D1R2[:,1] # DD-Paircount
norm=data_D1R2[0,2] # Normalization

for i, fn in enumerate(fns_D1R2[1:]):
    data_D1R2=np.loadtxt(fn) # Read in data_RR
    d1r2+=data_D1R2[:,1] # Add paircount
    norm+=data_D1R2[0,2] # Add normalization

d1r2/=norm


# Combine D2R1 if cross-correlations

if (not is_auto):
    data_D2R1=np.loadtxt(fns_D2R1[0])
    thetas=data_D2R1[:,0] # Bin centers for omega
    d2r1=data_D2R1[:,1] # DD-Paircount
    norm=data_D2R1[0,2] # Normalization

    for i, fn in enumerate(fns_D2R1[1:]):
        data_D2R1=np.loadtxt(fn) # Read in data_RR
        d2r1+=data_D2R1[:,1] # Add paircount
        norm+=data_D2R1[0,2] # Add normalization

    d2r1/=norm

# Calculate omega
if is_auto: #Auto correlations: Use Landy-Szalay
    omega = (dd-2*d1r2+rr)/rr
    omega=np.nan_to_num(omega) # Set to 0 if rr=0
else: #Cross correlations: Use Szalay-Szapudi
    omega = (dd-d1r2-d2r1+rr)/rr
    omega=np.nan_to_num(omega)


# Output
result=np.column_stack((thetas, omega))
np.savetxt(fn_out, result)

    
