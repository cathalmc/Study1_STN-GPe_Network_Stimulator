#!/bin/bash -l

# Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cores per node required
#SBATCH -n 10

# Set the walltime of the job to 10 mins (format is hh:mm:ss)
#SBATCH -t 6:00:00

# E-mail on begin (b), abort (a), and end (end) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=14369856@ucdconnect.ie

# Specify the jobname
#SBATCH --job-name=SWS

# Specify the error and output file names
#SBATCH --error="ErrorSW.out"
#SBATCH --output="OutputSW.out"

# Setup the environment
module load anaconda
module load gcc openmpi/3.1.4
conda activate --stack mpynn4

# Change to model working directory
model_dir="${HOME}/storage/ClusterRuns"

cd ${model_dir}

# Setup the model script
model_script=Commander.py
model="${model_dir}/${model_script}"

SLURM_NTASKS=6
NETWORK=

for (( i=0; i < ${SLURM_NTASKS} ; i++ )); do
  python ${model} ${SLURM_NTASKS} $i ${NETWORK} &
  if [ $(($i%9)) -eq 8 ]; then 
  sleep 1000
  fi

done
wait