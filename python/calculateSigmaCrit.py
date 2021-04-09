## @package calculateSigmaCrit.py
# Python script that calculates the critical surface mass density averaged
# over the source distribution
# Uses the Astropy.Cosmology package and Flat-LambdaCDM
# Source redshift distribution is read from file
#
# @param argv[1] Float, H_0 [km/s/Mpc]
# @param argv[2] Float, Omega_m
# @param argv[3] Filename of n(z) of sources
# @param argv[4] Float, lowest lens redshift
# @param argv[5] Float, highest lens redshift
# @param argv[6] Int, number of redshift bins
# @param argv[7] Outputfilename for \f$\bar{\Sigma}_\rm{crit}(z_l)$\f

import sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
from astropy.constants import c

## Calculates inverse sigma crit (without prefactor)
# @return \f $\frac{D_A(z_l, z_s) D_A(z_l)}{D_A(z_s)}$\f [Mpc]  
# @param cosmo Cosmology object (FlatLambdaCDM)
# @param z_l Lens redshift (float)
# @param z_s Source redshift (float)
def sigma_crit_inv(cosmo, z_l, z_s):
    D_s = cosmo.angular_diameter_distance(z_s)
    D_l = cosmo.angular_diameter_distance(z_l)
    D_ls = cosmo.angular_diameter_distance_z1z2(z_l, z_s)
    return (D_l*D_ls/D_s).to('Mpc')
    


# Parse cmd line
H0 = float(sys.argv[1])
Omm = float(sys.argv[2])
fn_nz = sys.argv[3]
z_low = float(sys.argv[4])
z_high = float(sys.argv[5])
z_nbins = int(sys.argv[6])
fn_out = sys.argv[7]

# Read in n(z)
data=np.loadtxt(fn_nz)
z_sources=data[:,0]
n_z=data[:,1]
delta_z_sources=data[1,0]-data[0,0]

# Create cosmology
cosmo=FlatLambdaCDM(H0=H0, Om0=Omm)

# Get constants
c_conv=c.to('Mpc/s') #speed-of-light in Mpc/c
G_conv=G.to('Mpc3/(Msun s2)')

# Create lens redshift bins
z_bins=np.linspace(z_low, z_high, z_nbins)
sig_crit=np.zeros(z_nbins)

# Calculate SigCrit
for i,z_l in enumerate(z_bins):
    for j,z_s in enumerate(z_sources):
        if(z_s>z_l):
            sig_crit[i]+=delta_z_sources*n_z[j]*sigma_crit_inv(cosmo, z_l, z_s).value


# Inverting
sig_crit=1/sig_crit
# Prefactor for Sigma_crit [Msun/Mpc]
prefactor=(c_conv*c_conv/4/np.pi/G_conv)
sig_crit*=prefactor.value

# Output
result=np.column_stack((z_bins, sig_crit))
np.savetxt(fn_out, result)
