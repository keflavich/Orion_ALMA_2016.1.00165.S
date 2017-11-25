#!/bin/sh

#This script is ment to be set in the COMMAND variable
#in the configure file to submit.  That submit script will create the
#clusterspec file for us in the WORK_DIR we specified in the configure file.

WORK_DIR='/lustre/aginsbur/orion/2016.1.00165.S/imaging'
cd ${WORK_DIR}
REDUCTION_DIR='/lustre/aginsbur/sgrb2/2016.1.00550.S/reduction'
REDUCTION_DIR=${WORK_DIR}

#export MPICASA=/home/casa/packages/RHEL6/release/casa-release-5.1.0-74/bin/mpicasa
#export CASA=/home/casa/packages/RHEL6/release/casa-release-5.1.0-74/bin/casa
#$MPICASA -n 4 -machinefile $PBS_NODE_FILE $CASA --nogui --nologger --log2term -c "scriptForImaging_lines_B3_TEST.py"

export CASAPATH=/home/casa/packages/RHEL6/release/casa-release-5.1.0-74
export PATH=${CASAPATH}/bin:$PATH
echo "PBS_NODEFILE = ",$PBS_NODEFILE
cat $PBS_NODEFILE
#mpicasa -machinefile $PBS_NODEFILE casa --nogui --nologger --log2term -c "${WORK_DIR}/scriptForImaging_lines_B3_TEST.py"
mpicasa -machinefile $PBS_NODEFILE casa --nogui --nologger --log2term -c "${WORK_DIR}/scriptForImaging_lines_B3.py"

#xvfb-run casapy --nogui -c "spwlist='0'; $SCRIPTPATH/fullcube_r0.py"
#xvfb-run casa --nogui -c "spwlist='0'; execfile('$SCRIPTPATH/fullcube_r0.py')"
