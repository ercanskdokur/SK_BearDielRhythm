#!/bin/bash
# ==============================================================================
# sub_R2_audit.sh — Audit-response batch: A1/A3/B1.3/B2/B4-B6/B8/C10/D1/D2/E2/E3
# ==============================================================================
set -uo pipefail
P=/path/to/project/BrownBearDielAct
SIF=$P/../skbear_r.sif

sub () {  # $1=script $2=name $3=walltime $4=mem $5=cpu $6=RESP(""=none) $7=extra_env
  local rtag="${6:-none}"
  sbatch --parsable \
    -J "$2_${rtag:0:3}" -p long -c "$5" --mem="$4" -t "$3" \
    -o "$P/logs/${2}_${rtag}_%j.out" -e "$P/logs/${2}_${rtag}_%j.err" \
    --wrap "unset LD_LIBRARY_PATH; module purge; module load singularity && ${6:+RESP=$6} TIME_AXIS=solar_double ${7:-} singularity exec --bind /path/to/project,/path/to/scratch $SIF Rscript $P/scripts/$1"
}

echo "== SINGLE-response jobs =="
echo "A3 moon(13)      $(sub 13_nocturnality_moon.R  a13  06:00:00 16G 4 '' '')"
echo "E3+B4 conflict   $(sub 15_conflict_windows.R   a15  02:00:00 16G 4 '' '')"
echo "D1 prior(09a)    $(sub 09a_prior_sensitivity.R a09a 10:00:00 48G 4 '' '')"
echo "B1.3 thresh(44)  $(sub 44_threshold_sens.R     a44  18:00:00 48G 4 '' '')"
echo "B8 blup(46)      $(sub 46_conflict_blup.R       a46  04:00:00 24G 4 '' '')"
echo "B2 loboSE(47)    $(sub 47_lobo_pairedSE.R       a47  02:00:00 32G 4 '' '')"

echo "== TWO-response jobs =="
for R in intensity duration; do
  echo "-- $R --"
  echo "C10 ppc(07b)    $(sub 07b_diagnostics_fixed.R c07b 05:00:00 48G 4 $R '')"
  echo "A1 collin(45)   $(sub 45_collinearity_sens.R  a45  12:00:00 32G 4 $R '')"
  echo "B4/5 tensor(42) $(sub 42_tensor_H3.R          a42  22:00:00 48G 4 $R '')"
  echo "B6 mundlak(43)  $(sub 43_mundlak.R            a43  14:00:00 32G 4 $R '')"
  echo "D2 beardayRI12b $(sub 12b_H3_bearday_RI.R     a12b 12:00:00 48G 4 $R '')"
  echo "D2 thinAR1 12c  $(sub 12c_H3_thinned_AR1.R    a12c 10:00:00 32G 4 $R '')"
  echo "E2 dump(29)     $(sub 29_dump_partitioning.R  a29  14:00:00 48G 4 $R '')"
  echo "E2 risk(30)     $(sub 30_conflict_risk_map.R  a30  10:00:00 48G 4 $R '')"
done
echo "== ALL JOBS SUBMITTED =="
