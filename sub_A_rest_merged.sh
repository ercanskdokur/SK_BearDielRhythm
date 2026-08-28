#!/bin/bash
# ==============================================================================
# sub_A_rest_merged.sh — REMAINDER of STAGE A as a SINGLE job (03 -> 03b -> 04 -> 05p)
# ==============================================================================
set -euo pipefail

P=/path/to/project/BrownBearDielAct
SIF=/path/to/project/skbear_r.sif

R_EXEC="singularity exec --bind /path/to/project,/path/to/scratch $SIF Rscript"

JID=$(sbatch --parsable \
  -J a_rest -p long -c 8 --mem=64G -t 48:00:00 \
  -o "$P/logs/a_rest_%j.out" -e "$P/logs/a_rest_%j.err" \
  --wrap "set -e; module load singularity; \
          echo \"NODE=\$(hostname) START=\$(date)\"; \
          $R_EXEC $P/scripts/03_covariate_extraction.R; \
          echo \"--- 03 DONE \$(date)\"; \
          $R_EXEC $P/scripts/03b_era5_hourly.R; \
          echo \"--- 03b DONE \$(date)\"; \
          $R_EXEC $P/scripts/04_decorrelate.R; \
          echo \"--- 04 DONE \$(date)\"; \
          $R_EXEC $P/scripts/05_prep_modeldata.R; \
          echo \"--- 05p DONE \$(date) — STAGE A COMPLETE\"")

echo "STAGE A (remainder) submitted as a SINGLE job: $JID"
echo "$JID" > "$P/logs/.last_stage_A"
echo
echo "Watch   : squeue -j $JID -o '%.10i %.9T %.10M %.20S %.20R'"
echo "Estimate: scontrol show job $JID | grep StartTime"
echo "Log     : tail -f $P/logs/a_rest_${JID}.out"
echo
echo "When done: bash sub_B_fits.sh   (NO argument — do NOT chain a dependency, so"
echo "           stage B starts accruing age the moment it is submitted; see point 1)"
