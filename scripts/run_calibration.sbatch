#!/bin/bash
#SBATCH --job-name=sojourn_calib
#SBATCH --output=logs/calib_%A_%a.out
#SBATCH --error=logs/calib_%A_%a.err
#SBATCH --array=0-9
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=PARTITION_NAME

# Ten array tasks x 300 replications = 3000 pooled replications.
# Pure NumPy, CPU only. Set --partition (and --account if your site
# requires one) to values appropriate for your cluster.

module load python 2>/dev/null || true
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs results

python src/sojourn_calib.py ${SLURM_ARRAY_TASK_ID} 300 \
    > results/calib_result_${SLURM_ARRAY_TASK_ID}.json

echo "task ${SLURM_ARRAY_TASK_ID} complete"
