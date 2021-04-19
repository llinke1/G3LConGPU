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

# Switch, true if autocorrelations, false if crosscorrelations
IS_AUTO=$4

# Switch, true if physical Gtilde is calculated
IS_PHYS=$5

# File with Tile Names
TILES=$6

# Number of tiles
NTILES=$7

# Redshiftweight width (if 0: No weighting)
SIGMA=$8

# File with Critical Surface Mass Density
FILE_SIGMACRIT=$9

# Maximal Number of Parallel Threads (Fit this to the number of processors)
MAX_JOBS=${10}

# Designator for galaxy files
NAME_OBJECTS=${11}

# Designator for mock files
NAME_SOURCES=${12}

# Switch if true Jackknifing is calculated
DO_JACKKNIFING=${13}

# Switch, if true paircount is calculated on GPU (brute-force), else on CPU (kdtrees)
GPU=${14}

# Folder with data files for first sample
DIR_DATA1=${15}

if [ "$IS_AUTO" = true ];
then   # Folder with data files for second sample
    DIR_DATA2=${16}
else
    DIR_DATA2=$DIR_DATA1
fi

# File with Critical Surface Mass Density
FILE_SIGMACRIT=${17}

# File with angular diameter distance
FILE_DA=${18}

# Binning for Gtilde
THETA_MIN=0.15
THETA_MAX=79.9 #199.75 #319.6 (MR) #79.9 (small tiles)
NBINS=128
R_MIN=0.02 #[Mpc]
R_MAX=40 #[Mpc] #70000000

# File with Angular Correlation Function
FILE_OMEGA=$DIR_PRODUCTS/omega/all.omega.dat


##################### Calculate Gtilde for Each Tile ##########################

if [ "$GPU" == 1 ];
then
    echo "Calculate Gtilde on GPU"
    BIN_GTILDE=$DIR_BIN/calculateGtilde_gpu.x
else
    echo "Calculate Gtilde on CPU (kd-Trees)"
    BIN_GTILDE=$DIR_BIN/calculateGtilde_kdtrees.x
fi

# Create folder for gtilde
mkdir -p $DIR_PRODUCTS/gtilde

# Do Gtilde calculation
for tile in $(awk 'NR>1 {print $1}' $TILES);
do
    if [ "$IS_PHYS" -gt 0 ]
    then
	echo ">GtildePhys for $tile | $(date)"
   
       $BIN_GTILDE $DIR_DATA1/$tile.$NAME_SOURCES.dat $DIR_DATA1/$tile.$NAME_OBJECTS.dat $DIR_DATA2/$tile.$NAME_OBJECTS.dat $FILE_OMEGA $THETA_MIN $THETA_MAX $R_MIN $R_MAX $NBINS $SIGMA $FILE_SIGMACRIT $FILE_DA $IS_PHYS  > $DIR_PRODUCTS/gtilde/$tile.gtilde_phys_single.dat
    else
	echo ">Gtilde for $tile | $(date)"
   
       $BIN_GTILDE $DIR_DATA1/$tile.$NAME_SOURCES.dat $DIR_DATA1/$tile.$NAME_OBJECTS.dat $DIR_DATA2/$tile.$NAME_OBJECTS.dat $FILE_OMEGA $THETA_MIN $THETA_MAX $R_MIN $R_MAX $NBINS $SIGMA $FILE_SIGMACRIT $FILE_DA $IS_PHYS  > $DIR_PRODUCTS/gtilde/$tile.gtilde_single.dat
    fi
done

###################### Combine Gtilde for all tiles #########################



echo ">Combining Gtilde | $(date)"
if [ "$IS_PHYS" -gt 0 ]
then
    python $DIR_PYTHON/combineGtilde.py 0 $NTILES $DIR_PRODUCTS/gtilde/all.gtilde_phys.dat $DIR_PRODUCTS/gtilde/*.gtilde_phys_single.dat 
else
  python $DIR_PYTHON/combineGtilde.py 0 $NTILES $DIR_PRODUCTS/gtilde/all.gtilde.dat $DIR_PRODUCTS/gtilde/*.gtilde_single.dat
fi

############### Combine Gtilde for each jackknife sample ####################

if [ "$DO_JACKKNIFING" -gt 0 ]
then
    echo ">Combining Gtilde for Jackknifing | $(date)"
    
    # Iterator Variable for Jackknife Samples
    I=1
    
    
    while [ $I -lt $NTILES ]; # Go through all Jackknifes
    do
	# Iterator Variable for Thread Number
	NUMBER_JOBS=0
	
	echo "Jackknifing for $I to $J+$MAX_JOBS | $(date)"
	
	while [ $NUMBER_JOBS -lt $MAX_JOBS ]; # Set Parallel Jobs
	do
	    if [ "$IS_PHYS" -gt 0 ]
	    then
		python $DIR_PYTHON/combineGtilde.py $I $NTILES $DIR_PRODUCTS/gtilde/jn_$I.gtilde_phys.dat $DIR_PRODUCTS/gtilde/*.gtilde_phys_single.dat  &
	    else
		python $DIR_PYTHON/combineGtilde.py $I $NTILES $DIR_PRODUCTS/gtilde/jn_$I.gtilde.dat $DIR_PRODUCTS/gtilde/*.gtilde_single.dat &
	    fi
	    ((NUMBER_JOBS++))
	    ((I++))
	    # Check if Number of Jackknife Samples is reached
	    if [ $I -ge $NTILES ];
	    then break
	    fi
	done
	wait # Wait for Jobs to finish
    done
fi

