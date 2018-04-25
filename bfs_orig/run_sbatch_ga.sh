#! /bin/bash

#SBATCH -p Lewis  # use the Lewis partition
#SBATCH -J run_ga  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00  # 2days -- twenty hour time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of cores (AKA tasks)
#SBATCH --mem=8GB # memory ~8GB
#SBATCH --exclusive # exlcusive node access to each job

## notifications
#SBATCH --mail-user=agtk4@mail.umsl.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#

# Commands here run only on the first core
echo "$(hostname), reporting for duty."

# Commands with srun will run on all cores in the allocation
./erbga $1 erbgaq-$1.out
