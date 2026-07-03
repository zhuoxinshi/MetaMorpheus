#!/usr/bin/env bash
# DDA the CONSENSUS-PAPER way: precursor feature tracing -> _ms1.feature -> mzLib FromFile precursor assembly.
# Consensus-traces the DDA MS1, writes a *_ms1.feature file, and assembles each MS2's precursor via
# ms2.GetIsolatedMassesAndCharges(ms1, FromFileDeconvolutionParameters) (the real MetaMorpheus join), then
# attaches the real fragments -> MGF -> td_pseudoMS2 search.
# (Uses the SearchDdaWithFromFileConsensusFeatures_FromEnv test harness.)
#
# Usage:  ./05_dda_fromfile_search.sh  <DDA.raw>  <label>
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"; source "$DIR/config.sh"
RAW="$1"; LABEL="$2"
fdir="$OUT/$LABEL"; mkdir -p "$fdir"
mz="$fdir/$(basename "${RAW%.raw}").mzML"
[ -f "$mz" ] || "$TRFP" -i "$RAW" -b "$mz"
mgf="$fdir/dda_fromfile.mgf"
DDAFF_MZML="$mz" DDAFF_OUT="$mgf" \
  dotnet test "$TESTPROJ" --no-build -c Debug --filter "FullyQualifiedName~SearchDdaWithFromFileConsensusFeatures_FromEnv"
sdir="$fdir/dda_fromfile_search"; mkdir -p "$sdir"
MM_LADDER_BONUS=0 dotnet "$CMD" -t "$TOMLS/td_pseudoMS2_search.toml" -d "$DB" -s "$mgf" -v minimal -o "$sdir"
echo "RESULT $LABEL DDA-FromFile-consensus (total ptm)= $(count_ids "$sdir")  spectra=$(grep -c 'BEGIN IONS' "$mgf")"
