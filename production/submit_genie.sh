#!/usr/bin/env bash
set -euo pipefail
# Script to simplify GENIE job submission to the grid

##### Parse any command-line options passed to this script ######
#
# Currently, the allowed options for invoking this script are
#
# -p, --prod:  Submit jobs with the VOMS Production role (rather than the
#              usual Analysis role). The job will be run on the grid
#              using the uboonepro account instead of your normal user
#              account. Job submission will fail unless you have permission
#              to submit production jobs.

# Number of expected command-line arguments
NUM_EXPECTED=8

# Process supported leading options in order until we see positional args or --.
while [ "$#" -gt 0 ]; do
    case "$1" in
        -p|--prod)
            USE_PRODUCTION="y"
            shift
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1"
            exit 2
            ;;
        *)
            break
            ;;
    esac
done

if [ "$#" -ne "$NUM_EXPECTED" ]; then
  echo "Usage: ./submit_genie.sh [--prod] NUM_JOBS JOB_NAME ARGS_FILE SPLINES_FILE FLUX_FILE FLUX_HIST_NAME OUTPUT_DIRECTORY GIT_CHECKOUT"
  exit 1
fi

NUM_JOBS=$1
STEM=$2
ARGS_FILE=$3
SPLINES_FILE=$4
FLUX_FILE=$5
FLUX_HIST_NAME=$6
OUTPUT_DIRECTORY=$7
GIT_CHECKOUT=$8

GRID_RESOURCES_DIR="${GRID_RESOURCES_DIR:-/pnfs/uboone/persistent/users/$(whoami)/grid}"

# Check if the USE_PRODUCTION variable is set
if [ -n "${USE_PRODUCTION+x}" ]; then
  echo "Production job will run as uboonepro"
  PROD_OPTION="--role=Production"
else
  echo "Analysis job will run as $(whoami)"
  PROD_OPTION="--role=Analysis"
fi

#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
#setup fife_utils

MY_GENIE_ARGS="$(<"$ARGS_FILE")"

jobsub_submit -G uboone --disk=60GB --expected-lifetime=8h -N "$NUM_JOBS" \
  $PROD_OPTION -f "$SPLINES_FILE" -f "$FLUX_FILE" \
  -e FLUX_HIST_NAME="${FLUX_HIST_NAME}" -d OUTPUT "$OUTPUT_DIRECTORY" \
  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE \
  --singularity-image \
    /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest \
  "file://${GRID_RESOURCES_DIR}/genie_grid.sh" "$STEM" "$(basename "$SPLINES_FILE")" \
  "$(basename "$FLUX_FILE")" "${GIT_CHECKOUT}" ${MY_GENIE_ARGS}
