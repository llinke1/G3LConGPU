#!/bin/bash

# Working Folder for Results, needs to already exist
DIR_PRODUCTS="../../test_triplecount/"

# Folder with C++ executables
DIR_BIN="../bin/"

# Files
FILE1="../../test_triplecount/galaxies1_very_large.dat"
FILE2="../../test_triplecount/galaxies2_very_large.dat"
FILE3="../../test_triplecount/galaxies3_very_large.dat"


# Binning
NBINS=15
R_MIN=0.02 #[Mpc]
R_MAX=2000 #[Mpc]

# File with Comoving distance
FILE_DCOM="comovingDist.dat"

mkdir -p $DIR_PRODUCTS

$DIR_BIN/calculateTriplecount_gpu.x $FILE1 $FILE2 $FILE3 $R_MIN $R_MAX $NBINS $FILE_DCOM > $DIR_PRODUCTS/ddd.dat

