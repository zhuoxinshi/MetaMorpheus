#!/usr/bin/env bash
# ISD: consensus tracing -> precursor-fragment grouping -> pseudoMS2 MGF -> top-down search -> ID count.
# Sweeps the precursor-fragment correlation cutoff. Optional: fragment mass-tracing only (no charge aggregation).
#
# Usage:  ./01_isd_consensus_search.sh  <ISD.raw or .mzML>  <label>  [corrList]  [fragAgg true|false]
# e.g.    ./01_isd_consensus_search.sh  $RAWDIR/09-18-25_YC_..._preFilter..._rep1.raw  YC  "0,0.4,0.6,0.7"  true
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"; source "$DIR/config.sh"
IN="$1"; LABEL="$2"; CORRS="${3:-0}"; FRAGAGG="${4:-true}"
fdir="$OUT/$LABEL"; mkdir -p "$fdir"

# (1) generate the pseudoMS2 MGF(s) via the env-driven test harness (builds consensus XICs once).
#     ISD_PREC_MINMASS/ISD_PREC_MINCHARGE filter the PRECURSOR channel to intact proteoforms (>=3 kDa, >=3 charges).
ISD_MZML="$IN" ISD_OUTDIR="$fdir" ISD_CORRS="$CORRS" ISD_FRAG_AGG="$FRAGAGG" \
  ISD_PREC_MINMASS=3000 ISD_PREC_MINCHARGE=3 \
  dotnet test "$TESTPROJ" --no-build -c Debug --filter "FullyQualifiedName~ConsensusIdSweep_FromEnv"

# (2) search each MGF with the top-down pseudoMS2 config and count IDs.
IFS=',' read -ra CS <<< "$CORRS"
for c in "${CS[@]}"; do
  mgf="$fdir/consensus_corr$c.mgf"; [ -f "$mgf" ] || { echo "$LABEL corr=$c NO_MGF"; continue; }
  sdir="$fdir/search_corr$c"; mkdir -p "$sdir"
  MM_LADDER_BONUS=0 dotnet "$CMD" -t "$TOMLS/td_pseudoMS2_search.toml" -d "$DB" -s "$mgf" -v minimal -o "$sdir"
  echo "RESULT $LABEL corr=$c  (total ptm)= $(count_ids "$sdir")  spectra=$(grep -c 'BEGIN IONS' "$mgf")"
done
