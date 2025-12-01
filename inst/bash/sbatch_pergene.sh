#!/bin/bash
#SBATCH --job-name=sba_TCGAhg19      # <-- change job name if desired
#SBATCH --account=ACF-EXAMPLE0001    # <-- change to your Slurm account
#SBATCH --output=sba_%j.out
#SBATCH --error=sba_%j.err
#SBATCH --time=24:00:00              # <-- change walltime for submission job
#SBATCH --mem=2G                     # <-- change memory for submission job
#SBATCH --cpus-per-task=1            # <-- change CPUs for submission job
#SBATCH --partition=campus           # <-- change to your partition
#SBATCH --qos=campus                 # <-- change to your QoS

# sbatch_pergene.sh
# Description:
#   Submit a higher-level job that, in turn, calls submit_pergene.sh
#   to schedule per-gene jobs.
#
# Usage:
#   sbatch sbatch_pergene.sh genelist_file study rscript_file submit_script
#
# Example:
#   sbatch sbatch_pergene.sh \
#     TCGAhg19_genelist.txt \
#     LUAD \
#     TCGAhg19_DR.R \
#     submit_pergene.sh


### ==== User-configurable settings (inside script body) ================== ###

R_MODULE="R/4.4.1"   # <-- change to your R module (or set empty if not using modules)

### ====================================================================== ###


if [ $# -ne 4 ]; then
  echo "Usage: sbatch $0 genelist_file study rscript_file submit_script"
  exit 1
fi

genelist_file="$1"
study="$2"
rscript_file="$3"
submit_script="$4"

if [ -n "${R_MODULE}" ]; then
  module load "${R_MODULE}"
fi

echo "Starting submission at: \$(date)"
echo "Working directory: \$(pwd)"

echo "Gene list file: ${genelist_file}"
echo "Study: ${study}"
echo "Rscript file: ${rscript_file}"
echo "Submit script: ${submit_script}"

bash "${submit_script}" "${genelist_file}" "${study}" "${rscript_file}"

echo "Finished submission at: \$(date)"
