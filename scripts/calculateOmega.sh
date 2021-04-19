#!/bin/bash

###############################################################################
# Skripts calculating Omega from Data
#
# First calculates DD, DR, RR for all tiles
# Second combines DD, DR, RR of all tiles
# Third, using the Landy-Szalay Estimator to get Omega for auto-correlations
# or Szalay-Szapudi Estimator for cross-correlations
#
# Author: Laila Linke, llinke@astro.uni-bonn.de
###############################################################################

# Working Folder for Results, needs to already exist
DIR_PRODUCTS=$1

# Folder with C++ executables
DIR_BIN=$2

# Folder with Python scripts
DIR_PYTHON=$3

# Switch, true if autocorrelations, false if crosscorrelations
IS_AUTO=$4

# File with tile names
TILES=$5

# Number of tiles
NTILES=$6

# Redshiftweight width (if 0: No weighting)
SIGMA=$7

# Designator for galaxy-files
NAME_OBJECTS=$8

# Designator for mock files
NAME_MOCKS=$9

# Switch, if true paircount is calculated on GPU (brute-force), else on CPU (kdtrees)
GPU=${10}

# Folder with data files for first sample
DIR_DATA1=${11}

if [ "$IS_AUTO" == 1 ];
then   # Folder with data files for second sample
    DIR_DATA2=$DIR_DATA1
else
    DIR_DATA2=${12}
fi

# Binning
THETA_MIN=0.013208 #[arcmin]
THETA_MAX=95.758 #212.13  #339.4 #95.758 #[arcmin]
N_BINS=${13}

##################### Calculate Paircounts for each tile ######################

if [ "$GPU" == 1 ];
then
    echo "Calculate Paircount on GPU"
    BIN_PAIRCOUNT=$DIR_BIN/calculatePaircount_gpu.x
else
    echo "Calculate Paircount on CPU (kd-Trees)"
    BIN_PAIRCOUNT=$DIR_BIN/calculatePaircount_kdtrees.x
fi
    
# Create Folder for omega
mkdir -p $DIR_PRODUCTS/omega

# Go through all tiles
for tile in $(awk 'NR>1 {print $1}' $TILES);
do
    echo ">Paircount for $tile | $(date)"
    # Calculate DD
    $BIN_PAIRCOUNT $DIR_DATA1/"$tile".$NAME_OBJECTS.dat $DIR_DATA2/"$tile".$NAME_OBJECTS.dat $SIGMA $THETA_MIN $THETA_MAX $N_BINS  > $DIR_PRODUCTS/omega/"$tile".DD.dat
   

    # Calculate D1R2
    $BIN_PAIRCOUNT $DIR_DATA1/"$tile".$NAME_OBJECTS.dat $DIR_DATA2/"$tile".$NAME_MOCKS.dat $SIGMA $THETA_MIN $THETA_MAX $N_BINS > $DIR_PRODUCTS/omega/"$tile".D1R2.dat 

    # Calculate RR
   $BIN_PAIRCOUNT $DIR_DATA1/"$tile".$NAME_MOCKS.dat $DIR_DATA2/"$tile".$NAME_MOCKS.dat $SIGMA $THETA_MIN $THETA_MAX $N_BINS > $DIR_PRODUCTS/omega/"$tile".RR.dat


    if [ "$IS_AUTO" == 0 ]; #Only for cross-correlations
    then
	# Calculate D2R1
	$BIN_PAIRCOUNT $DIR_DATA2/"$tile".$NAME_OBJECTS.dat $DIR_DATA1/"$tile".$NAME_MOCKS.dat $SIGMA $THETA_MIN $THETA_MAX $N_BINS > $DIR_PRODUCTS/omega/"$tile".D2R1.dat
    fi

done

##################### Combine Paircounts to Omega #############################

python $DIR_PYTHON/calculateOmega.py $IS_AUTO $NTILES $DIR_PRODUCTS/omega/all.omega.dat $DIR_PRODUCTS/omega/*.DD.dat $DIR_PRODUCTS/omega/*.RR.dat $DIR_PRODUCTS/omega/*.D1R2.dat $DIR_PRODUCTS/omega/*.D2R1.dat  
