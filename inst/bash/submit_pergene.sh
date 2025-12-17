#!/bin/bash

# submit_pergene.sh
# Description:
#   Submit per-gene jobs to Slurm for TCGA hg19 (or similar) data.
#   Each job runs an R script for a single gene.
#
# Usage:
#   ./submit_pergene.sh genelist.txt STUDY /path/to/script.R
#
# Example:
#   ./submit_pergene.sh TCGAhg19_genelist.txt LUAD TCGAhg19_DR.R


### ==== User-configurable settings (edit here) ============================ ###

# Slurm settings (edit to match your cluster and job size)
ACCOUNT="your_account"       # <-- change to your Slurm account
PARTITION="your_partition"   # <-- change to your partition name
QOS="your_qos"               # <-- change to your QoS (if used)

JOB_NAME_PREFIX="sub_"       # prefix for per-gene job names
TIME="3:00:00"               # walltime per gene job (e.g., 1:00:00, 3:00:00)
MEMORY="8G"                  # memory per gene job (e.g., 8G, 16G)
CPUS=1                       # CPUs per gene job

# Maximum number of concurrent "sub_" jobs for this user
MAX_JOBS=15                  # limit number of running/queued sub_ jobs

# Module settings (edit to match your environment)
MODULE_PURGE=true
MODULES_TO_LOAD="gnu openmpi R/4.4.1"  # e.g., "R/4.3.1" or similar

### ====================================================================== ###


if [ $# -ne 3 ]; then
  echo "Usage: $0 genelist_file study rscript_file"
  exit 1
fi

genelist_file="$1"
study="$2"
rscript_file="$3"

mkdir -p logs

while read gene; do
  # Check the number of SLURM jobs submitted by this user with name starting with JOB_NAME_PREFIX
  while [ "$(squeue -u "$USER" | grep -c "^${JOB_NAME_PREFIX}")" -ge "$MAX_JOBS" ]; do
    echo "Max jobs (${MAX_JOBS}) reached. Waiting 30 seconds..."
    sleep 30
  done

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME_PREFIX}${gene}
#SBATCH --account=${ACCOUNT}
#SBATCH --output=logs/${gene}_%j.out
#SBATCH --error=logs/${gene}_%j.err
#SBATCH --time=${TIME}
#SBATCH --mem=${MEMORY}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --partition=${PARTITION}
#SBATCH --qos=${QOS}

if [ "${MODULE_PURGE}" = "true" ]; then
  module purge
fi

if [ -n "${MODULES_TO_LOAD}" ]; then
  module load ${MODULES_TO_LOAD}
fi

echo "Running gene: ${gene} for Study: ${study}"
echo "Rscript: ${rscript_file}"
echo "Started at: \$(date)"

Rscript "${rscript_file}" "${study}" "${gene}"

echo "Finished at: \$(date)"
EOF

done < "${genelist_file}"
