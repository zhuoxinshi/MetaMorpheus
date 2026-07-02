#!/usr/bin/env bash
# DDA (normal): standard top-down search of the real MS2 run. Converts .raw -> .mzML first to avoid the
# interactive Thermo-license prompt that crashes a non-interactive CMD run.
#
# Usage:  ./02_dda_search.sh  <DDA.raw>  <label>
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"; source "$DIR/config.sh"
RAW="$1"; LABEL="$2"
fdir="$OUT/$LABEL"; mkdir -p "$fdir"
mz="$fdir/$(basename "${RAW%.raw}").mzML"
[ -f "$mz" ] || "$TRFP" -i "$RAW" -b "$mz"
sdir="$fdir/dda_search"; mkdir -p "$sdir"
dotnet "$CMD" -t "$TOMLS/dda_topdown_search.toml" -d "$DB" -s "$mz" -v minimal -o "$sdir"
echo "RESULT $LABEL DDA (total ptm)= $(count_ids "$sdir")"
