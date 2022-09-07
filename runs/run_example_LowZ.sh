#############################################################################
# Run for measurement of NNMap and NNMperp
#
# @author Laila Linke, llinke@astro.uni-bonn.de
#############################################################################


# Folder to which output shall be written
DIR_PRODUCTS="/vol/aibn238/data1/llinke/LowZ_IA/measurements/"
#"../../products_MR_m1m1_small_noZ/"

# Folder containing first lens and source files
DIR_DATA1="/vol/aibn238/data1/llinke/LowZ_IA/data/"#"../../MR/smallTiles_mass1/"

# Folder containing second lens files
DIR_DATA2="/vol/aibn238/data1/llinke/LowZ_IA/data/" #"../../MR/smallTiles_mass1/"


if [ $DIR_DATA1 == $DIR_DATA2 ]
then
    IS_AUTO=1
else
    IS_AUTO=0
fi

# Folder containing executables
DIR_BIN="../bin/"

# Folder containing Python scripts
DIR_PYTHON="../python"

# Folder containing Bash scripts
DIR_SCRIPTS="../scripts/"

# ASCII file with tile names
TILES=MR_small.tsb

# ASCII file with thetas for which NNMap shall be calculated [arcmin]
FILE_THETAS=thetas_intermed.dat

# ASCII file with Rs for which NNMap shall be calculated [Mpc]
FILE_RS=Rs_large.dat

SIGMA=0

# ASCII file containing Sigma Crit (source averaged)
FILE_SIGMACRIT="../../KV450/sigmaCrit/all.sigmaCrit.dat"

# ASCII file containing angular diameter distance
FILE_ANGULAR_DIST="../../KV450/D_A/angular_dist.dat"

# Prefix for lens galaxy files
NAME_OBJECTS=objects

# Prefix for random galaxy files
NAME_MOCKS=mocks

# Prefix for source galaxy files
NAME_SOURCES=sources

# Number of Jackknifes to be calculated
NUMBER_JN=1 # Match this to number of tiles

# Maximal number of parallel threads
MAX_JOBS=12 # Match this to number of cores

DO_JACKKNIFING=0

GPU=1 # Switch if GPU calculation

# ASCII file containing Sigma Crit (source averaged)
FILE_SIGMACRIT="../../KV450/sigmaCrit/all.sigmaCrit.dat"

# ASCII file containing angular diameter distance
FILE_DA="../../KV450/D_A/angular_dist.dat"

# ASCII file containing n_z
FILE_NZ="../../KV450/n_z/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB01t12_blindB_Nz.asc"

# Cosmology
H_0=73 # km/s/Mpc
OMEGA_M=0.25

# Binning
Z_MIN=0.001
Z_MAX=0.5
NBINS=128


# Create Outputfolder if it does not exist
mkdir -p $DIR_PRODUCTS


################ Calculate Omega ###############################################

# Calculate Omega
echo "Calculate Omega"
bash $DIR_SCRIPTS/calculateOmega.sh $DIR_PRODUCTS $DIR_BIN $DIR_PYTHON $IS_AUTO $TILES $NUMBER_JN $SIGMA $NAME_OBJECTS $NAME_MOCKS $GPU $DIR_DATA1 $DIR_DATA2 $NBINS

################ Calculate Sigma Crit ##########################################

#echo "Calculate SigmaCrit"
#python $DIR_PYTHON/calculateSigmaCrit.py $H_0 $OMEGA_M $FILE_NZ $Z_MIN $Z_MAX $NBINS $FILE_SIGMACRIT

################ Calculate D_A #################################################

#echo "Calculate D_A"
#python $DIR_PYTHON/calculateDA.py $H_0 $OMEGA_M $Z_MIN $Z_MAX $NBINS $FILE_ANGULAR_DIST

################ Calculate Gtilde (angular) ####################################

IS_PHYS=0

 Calculate Gtilde
bash $DIR_SCRIPTS/calculateGtildeLowZ.sh $DIR_PRODUCTS $DIR_BIN $DIR_PYTHON $IS_AUTO $IS_PHYS $TILES $NUMBER_JN $SIGMA $FILE_SIGMACRIT $MAX_JOBS $NAME_OBJECTS $NAME_SOURCES $DO_JACKKNIFING $GPU $DIR_DATA1 $DIR_DATA2 $FILE_SIGMACRIT $FILE_DA

################ Calculate NNMap (angular) #####################################

bash $DIR_SCRIPTS/calculateApertureStatistics.sh $DIR_PRODUCTS $DIR_BIN $DIR_PYTHON $IS_PHYS $FILE_RS $NUMBER_JN $MAX_JOBS $DO_JACKKNIFING


################ Calculate Gtilde (physical) ###################################

IS_PHYS=1

# Calculate Gtilde
#bash $DIR_SCRIPTS/calculateGtilde.sh $DIR_PRODUCTS $DIR_BIN $DIR_PYTHON $IS_AUTO $IS_PHYS $TILES $NUMBER_JN $SIGMA $FILE_SIGMACRIT $MAX_JOBS $NAME_OBJECTS $NAME_SOURCES $DO_JACKKNIFING $GPU $DIR_DATA1 $DIR_DATA2 $FILE_SIGMACRIT $FILE_DA

################ Calculate NNMap (physical) ####################################

#bash $DIR_SCRIPTS/calculateApertureStatistics.sh $DIR_PRODUCTS $DIR_BIN $DIR_PYTHON $IS_PHYS $FILE_RS $NUMBER_JN $MAX_JOBS $DO_JACKKNIFING


echo "Done"
