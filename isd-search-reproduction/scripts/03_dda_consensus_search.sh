#!/usr/bin/env bash
# DDA via CONSENSUS-TRACED precursors: consensus-trace the DDA MS1, pair each real DDA MS2 to its consensus
# precursor (isolation m/z + RT), keep the real fragments -> pseudoMS2 MGF -> same td_pseudoMS2 search.
# (Uses the SearchDdaWithConsensusPrecursors_FromEnv test harness.)
#
# Usage:  ./03_dda_consensus_search.sh  <DDA.raw>  <label>
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"; source "$DIR/config.sh"
RAW="$1"; LABEL="$2"
fdir="$OUT/$LABEL"; mkdir -p "$fdir"
mz="$fdir/$(basename "${RAW%.raw}").mzML"
[ -f "$mz" ] || "$TRFP" -i "$RAW" -b "$mz"
mgf="$fdir/dda_consensus.mgf"
DDA_MZML="$mz" DDA_OUT="$mgf" \
  dotnet test "$TESTPROJ" --no-build -c Debug --filter "FullyQualifiedName~SearchDdaWithConsensusPrecursors_FromEnv"
sdir="$fdir/dda_consensus_search"; mkdir -p "$sdir"
MM_LADDER_BONUS=0 dotnet "$CMD" -t "$TOMLS/td_pseudoMS2_search.toml" -d "$DB" -s "$mgf" -v minimal -o "$sdir"
echo "RESULT $LABEL DDA-consensus (total ptm)= $(count_ids "$sdir")  spectra=$(grep -c 'BEGIN IONS' "$mgf")"
