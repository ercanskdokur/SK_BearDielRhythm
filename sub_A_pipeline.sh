#!/bin/bash
# ==============================================================================
# sub_A_pipeline.sh — STAGE A: data pipeline (01 -> 02 -> 03 -> 03b -> 04 -> 05prep)
# ==============================================================================
# Usage:  bash sub_A_pipeline.sh
#
# NO MODEL FIT SHOULD BE SUBMITTED BEFORE THIS STAGE FINISHES.
# The chain runs from scratch into an empty project root so that no stale
# checkpoints / intermediate outputs (produced with the old, wrong time zone)
# are reused.
#
# The chain is SEQUENTIAL (--dependency=afterok): each step reads the previous
# step's output.
#   01 -> activity_clean / position_clean         (TZ FIX: UTC)
#   02 -> hourly_intensity (2 responses + SolarHour + t_period)
#   03 -> hourly_full (ctmm + rasters + E-OBS)
#   03b-> hourly_full + temp_hourly (downloads ERA5-Land!)
#   04 -> hourly_decorr + LandPC + VIF
#   05p-> model_data_50k.rds (shared 50K subsample)
# ==============================================================================
set -euo pipefail

P=/path/to/project/BrownBearDielAct
SIF=/path/to/project/skbear_r.sif

run_r () {  # $1=script  $2=job name  $3=walltime  $4=mem  $5=cpu  $6=dependency (may be empty)
  local dep=""
  [ -n "${6:-}" ] && dep="--dependency=afterok:$6"
  sbatch --parsable $dep \
    -J "$2" -p long -c "$5" --mem="$4" -t "$3" \
    -o "$P/logs/${2}_%j.out" -e "$P/logs/${2}_%j.err" \
    --wrap "module load singularity && singularity exec --bind /path/to/project,/path/to/scratch $SIF Rscript $P/scripts/$1"
}

echo "=== STAGE A: data pipeline ==="
J01=$(run_r 01_data_preparation.R    a01_prep    12:00:00 32G 4 "")            ; echo "01  -> $J01"
J02=$(run_r 02_response_variables.R  a02_resp    12:00:00 48G 4 "$J01")        ; echo "02  -> $J02"
J03=$(run_r 03_covariate_extraction.R a03_covs   24:00:00 64G 8 "$J02")        ; echo "03  -> $J03"
J3B=$(run_r 03b_era5_hourly.R        a03b_era5   12:00:00 32G 4 "$J03")        ; echo "03b -> $J3B  (DOWNLOADS ERA5)"
J04=$(run_r 04_decorrelate.R         a04_decorr  04:00:00 32G 4 "$J3B")        ; echo "04  -> $J04"
J05=$(run_r 05_prep_modeldata.R      a05_prep    04:00:00 32G 4 "$J04")        ; echo "05p -> $J05"

echo
echo "Chain submitted. Watch:  squeue -u \$USER"
echo "Last job: $J05  — once it finishes, sub_B_fits.sh can be run."
echo "$J05" > "$P/logs/.last_stage_A"
