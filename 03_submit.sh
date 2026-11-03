#!/bin/bash
# =========================================================================
# 03_submit.sh
#
# SLURM array job: each array task runs one parameter set of the design
# table (row <task id> of design.csv) via 02_run_batch.R.
#
# Setup (once, on the login node):
#   Rscript 00_make_design.R        # writes design.csv and design_check.png
#
# Smoke test (a few minutes):
#   DESIGN_FILE=design.csv N0_VAL=10000 GENS=200 N_STOCH_REPS=3 \
#     Rscript 02_run_batch.R 1
#   then inspect species_data_task1.csv and the printed peak clone count.
#
# Submit:
#   sbatch 03_submit.sh
#
# The design has 1024 rows. To run it in two batches of 512, submit with
# --array=1-512 and then --array=513-1024 (override on the command line,
# e.g. sbatch --array=513-1024 03_submit.sh). Keep every other setting
# identical between batches.
#
# Adjust --time, --mem-per-cpu and the array concurrency limit (%50) to
# your cluster, and set --mail-user if you want e-mail notifications.
# =========================================================================

#SBATCH --job-name=eco_evo
#SBATCH --array=1-1024%50
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.org
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/batch_%A_%a.out

mkdir -p logs results
module purge
module load R/4.3.3

# ---- simulation settings --------------------------------------------------
export NSPP=10; export GENS=1000; export RECORD_INTERVAL=25
export EARLY_WINDOW=50; export MAX_CLONES=2000000; export K_CEILING_MULT=100

# Initial abundance, the same for every species. Alternatively, set
# N0_FRAC (e.g. 0.5) to start each species at that fraction of its own
# single-species equilibrium r0 * K0; this overrides N0_VAL.
export N0_VAL=10000; export N0_FRAC=0

# Mutation effects: truncated normal, hard cap of +/-30% for r, K and alpha.
# The SD defaults to MAX_EFFECT / 3. A mean shift of 0 makes beneficial and
# deleterious mutations equally likely.
export MAX_EFFECT=0.3; export MUT_MEAN_SHIFT=0

# TRUE only for the null baseline (no mutations in the "with mutations" runs)
export NULL_TEST=FALSE

# Parameter set i uses row i of design.csv. Changing REPLICATES_PER_TASK or
# BASE_SEED changes which sets and random draws are used, so keep them
# the same across batches.
export REPLICATES_PER_TASK=1; export BASE_SEED=300000; export N_STOCH_REPS=100
export CHECKPOINT_INTERVAL=50
export MUTATION_FREQ_THRESHOLD=0.10

# Species-level draws within every set (log-uniform K and mu, uniform r).
# With a design file, K_MIN and K_MAX only fix the centre of the K range
# (geometric centre, ~3.2e5); its width in decades comes from design.csv.
export K_MIN=1e5; export K_MAX=1e6; export R_MIN=0.05; export R_MAX=0.2
export MU_MIN=1e-5; export MU_MAX=1e-3

# Community-level parameters (mean alpha, alpha CV, K spread) are read from here
export DESIGN_FILE=design.csv

# ---- run ------------------------------------------------------------------
LOG_FILE="logs/run_log_task${SLURM_ARRAY_TASK_ID}.txt"
if Rscript 02_run_batch.R "$SLURM_ARRAY_TASK_ID" > "$LOG_FILE" 2>&1; then
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully"
    mv "species_data_task${SLURM_ARRAY_TASK_ID}.csv" \
       "raw_freq_traj_task${SLURM_ARRAY_TASK_ID}.rds" \
       "mutation_log_task${SLURM_ARRAY_TASK_ID}.csv" \
       results/
else
    echo "ERROR: Task $SLURM_ARRAY_TASK_ID failed"; cat "$LOG_FILE"; exit 1
fi
