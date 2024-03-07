#!/bin/bash

###############################################################################
# Skripts calculating Gtilde from Data, Omega, Redshiftweight and Sigmacrit
#
# Author: Laila Linke, llinke@astro.uni-bonn.de
###############################################################################

# Working Folder for Results, needs to already exist
DIR_PRODUCTS=$1

# Folder with C++ executables
DIR_BIN=$2

# Folder with Python scripts
DIR_PYTHON=$3

# File with Tile Names
TILES=$4

# Maximal Number of Parallel Threads (Fit this to the number of processors)
MAX_JOBS=${5}

# Designator for galaxy files
NAME_OBJECTS=${6}

# Designator for mock files
NAME_SOURCES=${7}


# Switch, if true paircount is calculated on GPU (brute-force), else on CPU (kdtrees)
GPU=${8}

# Folder with data files for first sample
DIR_DATA1=${9}

DIR_DATA2=${10}

FLIP_E1=${11}
FLIP_E2=${12}

# Binning for Gplus / Gminus
THETA_MIN=0.02
THETA_MAX=79.9 #199.75 #319.6 (MR) #79.9 (small tiles)
NBINS=128


##################### Calculate Gplus/Gminus for Each Tile ##########################

if [ "$GPU" == 1 ];
then
    echo "Calculate Gplus / Gminus on GPU"
    BIN_GTILDE=$DIR_BIN/calculateGplusGminus_gpu.x
else
    echo "Calculate Gplus / Gminus on CPU (kd-Trees)"
    BIN_GTILDE=$DIR_BIN/calculateGplusGminus_kdtrees.x
fi

# Create folder for gtilde
mkdir -p $DIR_PRODUCTS/gplusgminus

# Do Gtilde calculation
for tile in $(awk 'NR>1 {print $1}' $TILES);
do
	echo ">Gtilde for $tile | $(date)"
   
    $BIN_GTILDE $DIR_DATA1/$tile.$NAME_OBJECTS.dat $DIR_DATA1/$tile.$NAME_SOURCES.dat $DIR_DATA2/$tile.$NAME_SOURCES.dat $THETA_MIN $THETA_MAX $NBINS $FLIP_E1 $FLIP_E2 > $DIR_PRODUCTS/gplusgminus/$tile.gplusgminus_single.dat
done

