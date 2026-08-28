#!/bin/bash
# ==============================================================================
# sub_B_fits.sh — STAGE B: 24 brms fits (6 hypotheses x 2 responses x 2 DUMPSPEC)
# ==============================================================================
# Usage:
#   bash sub_B_fits.sh              # if STAGE A is done (model_data_50k.rds exists)
#   bash sub_B_fits.sh 1234567      # submit with a dependency on a job
#
# ARRAY MAP (same order as HYP_GRID inside 05_fit_models.R):
#   expand.grid(hyp=c(H0,H1,H2,H3,H4,H5), resp=c(intensity,duration))
#   0=H0/int  1=H1/int  2=H2/int  3=H3/int  4=H4/int  5=H5/int
#   6=H0/dur  7=H1/dur  8=H2/dur  9=H3/dur 10=H4/dur 11=H5/dur
#
# WHY 24 FITS?
#   04_decorrelate finding: d2GarbageDump_km_sc has VIF=12.31, R2_on_others=0.9188
#   -> the study's focal variable is 92% explained by the other covariates. The
#   dominant shared axis is LandPC1 ("distance from everything"); the dump lies a
#   few km from town. Decision: fit BOTH versions and compare by LOO (sensitivity).
#     DUMPSPEC=raw   -> raw                 (unsuffixed file; PRIMARY)
#     DUMPSPEC=resid -> residualised on LandPC1+LandPC2 (_dumpresid suffix)
#
#   PRE-VALIDATED (same code as 05_fit_models.R):
#     VIF 12.32 -> 2.90 | 49,979 rows preserved | no NA | sd=1
#     correlation with raw = 0.485 -> genuinely different signal, not cosmetic
#     (76.5% of dump variance was explained by LandPC1+PC2)
#
# CRITICAL: both versions share the same response and the same ROWS; only the
#   CONTENT of the dump column changes -> elpd is DIRECTLY comparable.
#
# NOTE: H3 and H5 are heavy (4 by-smooths / random slopes) -> all get the same
# generous resources; light models finish early, which is fine.
# ==============================================================================
set -euo pipefail

P=/path/to/project/BrownBearDielAct
SIF=/path/to/project/skbear_r.sif
DEP=""
[ -n "${1:-}" ] && DEP="--dependency=afterok:$1"

if [ ! -f "$P/models/model_data_50k.rds" ] && [ -z "${1:-}" ]; then
  echo "ERROR: model_data_50k.rds missing. Complete STAGE A first." >&2
  exit 1
fi

submit_variant () {   # $1 = DUMPSPEC (raw | resid)
  sbatch --parsable $DEP \
    -J "b_$1" -p long -c 4 --mem=48G -t 72:00:00 \
    --array=0-11%6 \
    -o "$P/logs/b_${1}_%A_%a.out" -e "$P/logs/b_${1}_%A_%a.err" \
    --wrap "module load singularity && \
            TIME_AXIS=solar_double DUMPSPEC=$1 \
            singularity exec --bind /path/to/project,/path/to/scratch $SIF Rscript $P/scripts/05_fit_models.R"
}

J_RAW=$(submit_variant raw)
J_RES=$(submit_variant resid)

echo "STAGE B submitted — 24 fits, 2 arrays (each 0-11, at most 6 concurrent):"
echo "  DUMPSPEC=raw   -> $J_RAW   outputs: H*_<resp>.rds            (PRIMARY)"
echo "  DUMPSPEC=resid -> $J_RES   outputs: H*_<resp>_dumpresid.rds  (sensitivity)"
printf '%s\n%s\n' "$J_RAW" "$J_RES" > "$P/logs/.last_stage_B"
echo
echo "Watch: squeue -u \$USER"
echo "Log  : ls -t $P/logs/b_raw_* $P/logs/b_resid_* | head"
echo "Out  : ls $P/models/H*_intensity*.rds $P/models/H*_duration*.rds 2>/dev/null"
echo
echo "When done: bash sub_C_derived.sh  (derived analyses on the PRIMARY models)"
echo "           + raw vs resid LOO comparison -> the data decides which spec is reported"
