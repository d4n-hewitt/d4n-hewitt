#!/bin/bash -l

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request one hour of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=48:00:00

# 3. Request 1 gigabyte of RAM per core
#$ -l mem=1G

# 4. Set the name of the job.
#$ -N MFI_MdcTST

# 6. Select the MPI parallel environment and 12 processors.
#$ -pe smp 1

#$ -wd /home/zccadhe/Scratch/Xylene_MD/MFI/dcTST/Meta
# 7. Set the working directory to somewhere in your scratch space.  
# Replace <your userid> with your UCL userid.
~/RASPA_INSTALL/bin/simulate simulation.input
