#!/usr/bin/env bash
# GPTMD -> Search (2-task chain) for either an ISD pseudoMS2 MGF or a DDA .mzML.
# The GPTMD task augments the database with discovered mods; the SearchTask then searches the augmented DB.
#
# Usage:  ./04_gptmd.sh  isd  <ISD_pseudoMS2.mgf>   <label>
#         ./04_gptmd.sh  dda  <DDA.mzML>            <label>
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"; source "$DIR/config.sh"
MODE="$1"; SPEC="$2"; LABEL="$3"
fdir="$OUT/$LABEL"; mkdir -p "$fdir"
if [ "$MODE" = "isd" ]; then
  gptmd="$TOMLS/gptmd_isd.toml"; search="$TOMLS/td_pseudoMS2_search.toml"; sdir="$fdir/isd_gptmd"
else
  gptmd="$TOMLS/gptmd_dda.toml"; search="$TOMLS/dda_topdown_search.toml"; sdir="$fdir/dda_gptmd"
fi
mkdir -p "$sdir"
# multiple -t tomls run in sequence; the search consumes the GPTMD-augmented DB
MM_LADDER_BONUS=0 dotnet "$CMD" -t "$gptmd" "$search" -d "$DB" -s "$SPEC" -v minimal -o "$sdir"
echo "RESULT $LABEL $MODE+GPTMD (total ptm)= $(count_ids "$sdir")"
