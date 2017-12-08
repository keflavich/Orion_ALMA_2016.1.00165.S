#!/bin/sh

#This script is ment to be set in the COMMAND variable
#in the configure file to submit.  That submit script will create the
#clusterspec file for us in the WORK_DIR we specified in the configure file.

WORK_DIR='/lustre/aginsbur/orion/2016.1.00165.S/imaging'
cd ${WORK_DIR}
REDUCTION_DIR=${WORK_DIR}

export CASAPATH=/home/casa/packages/RHEL6/release/casa-release-5.1.0-74
export PATH=${CASAPATH}/bin:$PATH
echo "PBS_NODEFILE = ",$PBS_NODEFILE
cat $PBS_NODEFILE
mpicasa -machinefile $PBS_NODEFILE casa --nogui --nologger --log2term -c "${WORK_DIR}/test480to640_scriptForImaging_lines_B6.py"
