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

# File whith thetas for which N2Map shall be calculated
FILE_THETAS=$4


# Maximal Number of Parallel Threads (Fit this to the number of processors)
MAX_JOBS=$5

DO_JACKKNIFING=$6

TILES=$7


###############################################################################

echo ">Aperture Statistics w/o tesselation | $(date)"

for tile in $(awk 'NR>1 {print $1}' $TILES);
do
    echo Processing $tile
    $DIR_BIN/calculateNMM.x $DIR_PRODUCTS/gplusgminus/$tile.gplusgminus_single.dat $FILE_THETAS > $DIR_PRODUCTS/NNMap/$tile.NMM.dat

done