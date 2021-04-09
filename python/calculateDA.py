## @package calculateDA.py
# Python script that calculates the angular diameter distances [Mpc]
# Uses the Astropy.Cosmology package and Flat-LambdaCDM
#
# @param argv[1] Float, \f$H_0$\f [km/s/Mpc]
# @param argv[2] Float, \f$\Omega_m$\f
# @param argv[3] Float, lowest lens redshift
# @param argv[4] Float, highest lens redshift
# @param argv[5] Int, number of redshift bins
# @param argv[6] Outputfilename for \f$\bar{\Sigma}_\rm{crit}(z_l)$\f

import sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM

# Parse cmd line
H0 = float(sys.argv[1])
Omm = float(sys.argv[2])
z_low = float(sys.argv[3])
z_high = float(sys.argv[4])
z_nbins = int(sys.argv[5])
fn_out = sys.argv[6]

# Create cosmology
cosmo=FlatLambdaCDM(H0=H0, Om0=Omm)

# Create lens redshift bins
z_bins = np.linspace(z_low, z_high, z_nbins)

# Calculate angular diameter distance
D_A = cosmo.angular_diameter_distance(z_bins).value

# Output
result = np.column_stack((z_bins, D_A))
np.savetxt(fn_out, result)
