#! /bin/bash
echo hello world
### combine PBS standard output and error files
#PBS -j oe
#PBS -o $PBS_JOBNAME.o$PBS_JOBID

### change to working directory where you submit job 
cd $PBS_O_WORKDIR

### use PBS array id to specify chromosome
if [ "$PBS_ARRAYID" == "undefined" ] || [ "$PBS_ARRAYID" == "" ]; then
    TASK=""
else
    TASK="--chromosome $PBS_ARRAYID"
fi

R -q --vanilla --args ${args} $TASK < ${R}
