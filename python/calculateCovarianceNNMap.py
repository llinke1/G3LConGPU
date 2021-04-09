## @package calculateCovarianceNNMap.py
# Calculate Covariancematrix and Stddev of <N2Map> from Jackknife samples
#
# @param argv[1] Outputfilename for covariance matrix
# @param argv[2] Outputfilename for standard deviation
# @param argv[3:N] Jackknifesamples
# @author Laila Linke, llinke@astro.uni-bonn.de

import numpy as np
import sys

message="""
# calculateCovarianceNNMap.py: WRONG NUMBER OF ARGUMENTS
# Usage: calculateCovarianceNNMap.py Outputfile_Covariance_Matrix
#	 Outputfile_StandardDeviations NNMap_Jackknifes
#
# Example: calculateCovarianceNNMap 
#	   ../../products/NNMap/all.e_mode_covariance.dat
#	   ../../products/NNMap/all.e_mode_stddev.dat
#          ../../products/NNMap/jn*.N2Map.dat
"""

if(len(sys.argv) < 4): # Check Commandline
    print (message)
    exit()



# Reading in Filenames
filename_covariance=sys.argv[1] #Outputfile for Covariance Matrix
filename_variance=sys.argv[2] # Outputfile for Standarddeviations
filelist_jn=sys.argv[3:] # List of Jackknife files

number_jn=len(filelist_jn) # Number of Jackknifes

# Initialization
jackknife_1=np.loadtxt(filelist_jn[0]) # Read in first Jackknife Sample

average_jn=jackknife_1[:,3] # Average of all jackknifes
thetas=jackknife_1[:,0] # Theta1 of N2Map

# Calculate Average of Jackknifes
for f in filelist_jn[1:]:
        jackknife=np.loadtxt(f)
        average_jn+=jackknife[:,3]
        
average_jn/=number_jn

# Calculate Covariance Matrix

number_data=len(average_jn)

       
covariance=np.zeros((number_data, number_data))

for f in filelist_jn:
        data=np.loadtxt(f)
        jackknife=data[:,3]
        for i in range(number_data):
                for j in range(number_data):
                        covariance[i][j]+=(average_jn[i]-jackknife[i])*(average_jn[j]-jackknife[j])
                        
covariance*=(number_jn-1.0)/number_jn

# Save Covariance Matrix
np.savetxt(filename_covariance, covariance)

# Calculate Standarddeviation (Sqrt of Diagonal of Covariance Matrix)
variance=np.sqrt(np.diagonal(covariance))
a=np.column_stack((thetas,variance))

# Save Standarddeviation
np.savetxt(filename_variance,a)
