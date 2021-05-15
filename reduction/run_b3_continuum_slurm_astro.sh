#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --mem=300gb                     # Job memory request
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept
pwd; hostname; date

module load git

which python
which git

git --version
echo $?

cd /orange/adamginsburg/orion/2016.1.00165.S/imaging

scriptpath=/orange/adamginsburg/orion/Orion_ALMA_2016.1.00165.S/reduction

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
#export CASA=/orange/adamginsburg/casa/casa-release-5.6.0-60.el7/bin/casa

xvfb-run -d ${CASA}  --nogui --nologger -c "execfile('${scriptpath}/continuum_imaging_b3_sep17_2020.py')"
