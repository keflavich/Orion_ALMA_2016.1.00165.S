#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail	
#SBATCH --ntasks=8
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --nodes=1
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
pwd; hostname; date

module load git

which python
which git

git --version
echo $?

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
#export CASA=/orange/adamginsburg/casa/casa-release-5.6.0-60.el7/bin/casa

xvfb-run -d ${CASA}  --nogui --nologger -c "execfile('continuum_imaging_b6_sep17_2020.py')"
