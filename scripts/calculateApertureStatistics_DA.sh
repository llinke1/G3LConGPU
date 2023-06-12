#!/bin/bash

###############################################################################
# Skripts calculating <N2Map> from Gtilde including Tesselation
#
# Author: Laila Linke, llinke@astro.uni-bonn.de
###############################################################################

# Working Folder for Results, needs to already exist
DIR_PRODUCTS=$1

# Folder with C++ executables
DIR_BIN=$2

# Folder with Python scripts
DIR_PYTHON=$3

# Switch, true if physical Gtilde is calculated
IS_PHYS=$4

# File whith thetas for which N2Map shall be calculated
FILE_THETAS=$5

# Number of Jackknife Samples
NUMBER_JN=$6

# Maximal Number of Parallel Threads (Fit this to the number of processors)
MAX_JOBS=$7

DO_JACKKNIFING=$8

TILES=$9

#if [ "$IS_PHYS" -gt 0 ]
#then
#    GTILDE=gtilde_phys
#    N2MAP=N2Map_phys
#else
GTILDE=gtilde
N2MAP=N2Map_DA
#fi

###############################################################################

############ Function for Computing the Tesselation ###########################
# $1: all or jn_NUMBER
###############################################################################

function tesselation
{

    
    # Find all filled Bins and write to temporary file
    awk '($9+0) > 0 {print NR, $1, $2, $3}' $DIR_PRODUCTS/gtilde/$1.$GTILDE.dat > $DIR_PRODUCTS/gtilde/$GTILDE.filledBins_$1

	voro++ 0.1 100 0.1 100 0 6.2832 $DIR_PRODUCTS/gtilde/gtilde.filledBins_$1


#    if [ "$IS_PHYS" -gt 0 ]
#    then
#	voro++ 0.02 30 0.02 30 0 6.2832 $DIR_PRODUCTS/gtilde/gtilde_phys.filledBins_$1
#    else
	# Do tesselation and write volumes to temporary file.vol
#	voro++ 0.05 79.9 0.05 79.9 0 6.2832 $DIR_PRODUCTS/gtilde/gtilde.filledBins_$1
#    fi

    

    # Combine bin, gtilde and tesselation volumes in file
    
     awk '{print NR, $0}' $DIR_PRODUCTS/gtilde/$1.$GTILDE.dat > $DIR_PRODUCTS/gtilde/$1_numbered.$GTILDE.dat

     join -1 1 -2 1 <(sort -k1b,1 -n $DIR_PRODUCTS/gtilde/$1_numbered.$GTILDE.dat) <(sort -k1b,1 -n $DIR_PRODUCTS/gtilde/$GTILDE.filledBins_$1.vol)|  awk '{print $2, $3, $4, $14, $8, $9, $10}' > $DIR_PRODUCTS/gtilde/$1.$GTILDE.tesselated.dat

     #rm $DIR_PRODUCTS/gtilde/$GTILDE.filledBins_$1
     #rm $DIR_PRODUCTS/gtilde/$GTILDE.filledBins_$1.vol
     #rm $DIR_PRODUCTS/gtilde/$1_numbered.$GTILDE.dat

}

mkdir -p $DIR_PRODUCTS/NNMap

###############################################################################

echo ">Aperture Statistics w/o tesselation | $(date)"
echo Tiles: $TILES
echo Thetates: $FILE_THETAS
# for tile in $(awk 'NR>1 {print $1}' $TILES);
# do
#     echo Processing $tile
#     $DIR_BIN/calculateApertureStatistics.x $DIR_PRODUCTS/gtilde/$tile.gtilde_single.dat $FILE_THETAS 0 > $DIR_PRODUCTS/NNMap/$tile.$N2MAP.dat

# done

$DIR_BIN/calculateApertureStatistics.x $DIR_PRODUCTS/gtilde/all.gtilde_DA.dat $FILE_THETAS 0 > $DIR_PRODUCTS/NNMap/all.$N2MAP.dat


# if [ $DO_JACKKNIFING -gt 0 ];
# then
    ############### Calculate <N2Map> for each Jackknife Sample ##################

    echo ">Calculating Jackknifing | $(date)"
    
    # Iterator Variable for Jackknife Samples
    J=1
    
    while [ $J -lt $NUMBER_JN ]; # Go through all Jackknifes
    do
	# Iterator Variable for Thread Number
	NUMBER_JOBS=0
	
	echo ">Jackknifing for $J to $J+$MAX_JOBS | $(date)"
	
	    while [ $NUMBER_JOBS -lt $MAX_JOBS ]; # Set Parallel Jobs
	    do
	        $DIR_BIN/calculateApertureStatistics.x $DIR_PRODUCTS/gtilde/jn_$J.gtilde_DA.dat $FILE_THETAS 0 > $DIR_PRODUCTS/NNMap/jn_$J.$N2MAP.dat &
	        ((NUMBER_JOBS++))
	        ((J++))
	        # Check if Number of Jackknife Samples is reached
	    if [ $J -ge $NUMBER_JN ];
	    then break
	    fi
        done
	done
# fi

# ############ Do Tesselation for All Gtilde ####################################

#echo ">Tesselation | $(date)"

#tesselation all

# ################ Calculate <N2Map> for whole Sample ##########################

# mkdir -p $DIR_PRODUCTS/NNMap



#echo ">Aperture Statistics with tesselation | $(date)"
#$DIR_BIN/calculateApertureStatistics.x $DIR_PRODUCTS/gtilde/all.$GTILDE.tesselated.dat $FILE_THETAS 1 > $DIR_PRODUCTS/NNMap/all.$N2MAP.tesselated.dat

# ############ Do Tesselation for Jackknife Samples #############################


# if [ $DO_JACKKNIFING -gt 0 ];
# then
    
#     echo ">Tesselation for Jackknifing | $(date)"

#     # Iterator Variable for Jackknife Samples
#     I=1
    
#     while [ $I -lt $NUMBER_JN ]; # Go through all Jackknifes
#     do
# 	# Iterator Variable for Thread Number
# 	NUMBER_JOBS=0
	
# 	echo ">Started Tesselation for $I to $I+$MAX_JOBS | $(date)"
	
# 	while [ $NUMBER_JOBS -lt $MAX_JOBS ]; # Set Parallel Jobs
# 	do
# 	    tesselation jn_$I &
# 	    ((NUMBER_JOBS++))
# 	    ((I++))
# 	    # Check if Number of Jackknife Samples is reached
# 	    if [ $I -ge $NUMBER_JN ];
# 	    then break
# 	    fi
# 	done
# 	wait # Wait for Jobs to finish
#     done
    
#     ############### Calculate <N2Map> for each Jackknife Sample ##################

#     echo ">Calculating Jackknifing | $(date)"
    
#     # Iterator Variable for Jackknife Samples
#     J=1
    
#     while [ $J -lt $NUMBER_JN ]; # Go through all Jackknifes
#     do
# 	# Iterator Variable for Thread Number
# 	NUMBER_JOBS=0
	
# 	echo ">Jackknifing for $J to $J+$MAX_JOBS | $(date)"
	
# 	while [ $NUMBER_JOBS -lt $MAX_JOBS ]; # Set Parallel Jobs
# 	do
# 	    $DIR_BIN/calculateApertureStatistics.x $DIR_PRODUCTS/gtilde/jn_$J.$GTILDE.tesselated.dat $FILE_THETAS > $DIR_PRODUCTS/NNMap/jn_$J.$N2MAP.tesselated.dat &
# 	    ((NUMBER_JOBS++))
# 	    ((J++))
# 	    # Check if Number of Jackknife Samples is reached
# 	if [ $J -ge $NUMBER_JN ];
# 	then break
# 	fi
# 	done
# 	wait # Wait for Jobs to finish
#     done
    
# ############### Calculate Variance from Jackknife Samples #####################
    
#     echo ">Calculating Variance | $(date)"
    
#     python $DIR_PYTHON/calculateCovarianceNNMap.py $DIR_PRODUCTS/NNMap/all.e_mode_covariance_tesselated.dat $DIR_PRODUCTS/NNMap/all.e_mode_stddev_tesselated.dat $DIR_PRODUCTS/NNMap/jn*.$N2MAP.tesselated.dat

#     python $DIR_PYTHON/calculateCovarianceNNMperp.py $DIR_PRODUCTS/NNMap/all.b_mode_covariance_tesselated.dat $DIR_PRODUCTS/NNMap/all.b_mode_stddev_tesselated.dat $DIR_PRODUCTS/NNMap/jn*.$N2MAP.tesselated.dat
    
#     paste -d " " <(awk '{print $0}' $DIR_PRODUCTS/NNMap/all.N2Map_tesselated.dat ) <(awk '{print $2}' $DIR_PRODUCTS/NNMap/all.e_mode_stddev_tesselated.dat ) <(awk '{print $2}' $DIR_PRODUCTS/NNMap/all.b_mode_stddev_tesselated.dat ) >  $DIR_PRODUCTS/NNMap/all.$N2MAP.tesselated_variance.dat
    
# fi
