#!/bin/bash
# ==============================================================================
# sub_C_derived.sh — STAGE C: diagnostics, comparison, figures, derived analyses
# ==============================================================================
# Usage:
#   bash sub_C_derived.sh            # if STAGE B is done
#   bash sub_C_derived.sh 1234567    # depend on B's array id
#
# Each script runs separately for the two responses (RESP=intensity|duration);
# outputs are suffixed via tag_out() so they do NOT overwrite each other.
#
# 09c (grouped bear-CV) is the heaviest job: submit it separately, on a long queue.
# ==============================================================================
set -euo pipefail

P=/path/to/project/BrownBearDielAct
SIF=/path/to/project/skbear_r.sif
DEP=""
[ -n "${1:-}" ] && DEP="--dependency=afterok:$1"

sub () {  # $1=script $2=name $3=walltime $4=mem $5=cpu $6=RESP $7=extra_env
  sbatch --parsable $DEP \
    -J "$2_${6:0:3}" -p long -c "$5" --mem="$4" -t "$3" \
    -o "$P/logs/${2}_${6}_%j.out" -e "$P/logs/${2}_${6}_%j.err" \
    --wrap "module load singularity && RESP=$6 TIME_AXIS=solar_double ${7:-} \
            singularity exec --bind /path/to/project,/path/to/scratch $SIF Rscript $P/scripts/$1"
}

echo "=== STAGE C: derived analyses ==="
for R in intensity duration; do
  echo "--- RESP=$R ---"
  echo "  06b LOO         -> $(sub 06b_loo_extended.R        c06b  12:00:00 64G 4 $R)"
  echo "  06c dumpspec    -> $(sub 06c_dumpspec_compare.R    c06c  12:00:00 64G 4 $R)"
  echo "  07b diag+VPC    -> $(sub 07b_diagnostics_fixed.R   c07b  12:00:00 64G 4 $R)"
  echo "  10  viz M1-M5   -> $(sub 10_visualization.R        c10   04:00:00 32G 4 $R)"
  echo "  10b viz H3/H4   -> $(sub 10b_visualization_newmodels.R c10b 04:00:00 32G 4 $R)"
  echo "  14  individual  -> $(sub 14_individual_specialization.R c14 04:00:00 32G 4 $R)"
  echo "  32  level/time  -> $(sub 32_level_vs_timing.R      c32   08:00:00 48G 4 $R)"
  echo "  29  dump partit.-> $(sub 29_dump_partitioning.R    c29   24:00:00 48G 4 $R)"
  echo "  30  risk map    -> $(sub 30_conflict_risk_map.R    c30   12:00:00 64G 4 $R)"
  echo "  12a autocorr    -> $(sub 12a_autocorr_diag.R       c12a  04:00:00 32G 4 $R)"
done

# Response-independent (builds its own model)
echo "--- response-independent ---"
echo "  11  nocturnality (strict) -> $(sub 11_nocturnality_index.R c11s 12:00:00 48G 4 intensity 'NIGHT_DEF=strict')"
echo "  11  nocturnality (wide)   -> $(sub 11_nocturnality_index.R c11w 12:00:00 48G 4 intensity 'NIGHT_DEF=wide')"
echo "  15  conflict windows      -> $(sub 15_conflict_windows.R   c15  04:00:00 32G 4 intensity)"

echo
echo "SUBMIT SEPARATELY (very heavy, own queue):"
echo "  09c grouped bear-CV : bash sub_D_lobo.sh"
echo "  31  axis comparison : needs an H3 fit with TIME_AXIS=clock first (see header of 31)"
echo "  13  lunar control   : after 11 finishes"
echo "  16  phenotypes      : after 11 + 14 finish"
