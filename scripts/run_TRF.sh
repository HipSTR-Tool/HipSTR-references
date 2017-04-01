# Check the number of command line arguments
if [ $# != 3 ]
then
    echo "ERROR: This script requires three arguments: the FASTA file, output directory and the minimum score. Exiting..."
    exit 1
fi

# Input parameters
INPUT_FILE=$1
OUTPUT_DIR=$2
MIN_SCORE=$3

# Default parameters to use with Tandem Repeat Finder
MATCH_WT=2
MISMATCH_PEN=7
INDEL_PEN=7
P_MATCH=80
P_INDEL=10
MAX_PERIOD=500

NAME=$(basename $INPUT_FILE)

# Create output directory if it doesn't already exist
if [ ! -d $OUTPUT_DIR ]
then
    echo "Creating output directory $OUTPUT_DIR"
    mkdir $OUTPUT_DIR
fi

./trf409.legacylinux64 $INPUT_FILE $MATCH_WT $MISMATCH_PEN $INDEL_PEN $P_MATCH $P_INDEL $MIN_SCORE $MAX_PERIOD -h -d -l 6 -ngs > $OUTPUT_DIR/$NAME



