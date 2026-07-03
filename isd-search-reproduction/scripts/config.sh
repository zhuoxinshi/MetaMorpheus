#!/usr/bin/env bash
# Edit these paths for your machine, then `source config.sh` from the other scripts.
# All scripts assume you have BUILT this branch (Release CMD + Debug Test):
#   dotnet build MetaMorpheus/MetaMorpheus/CMD/CMD.csproj   -c Release
#   dotnet build MetaMorpheus/MetaMorpheus/Test/Test.csproj -c Debug

# Repo root of THIS MetaMorpheus checkout (contains MetaMorpheus/MetaMorpheus/...)
REPO="E:/CodeReview/ISD/MetaMorpheus"

CMD="$REPO/MetaMorpheus/CMD/bin/Release/net8.0/CMD.dll"
TESTPROJ="$REPO/MetaMorpheus/Test/Test.csproj"
TOMLS="$REPO/isd-search-reproduction/tomls"

# Yeast UniProt (reviewed) database used for the searches
DB="E:/CodeReview/ISD/Test_runs/09-17-25/uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml"

# Directory holding the .raw / .mzML acquisitions
RAWDIR="E:/CodeReview/ISD/Test_runs/09-17-25"

# Where to write generated MGFs + search outputs
OUT="E:/CodeReview/ISD/deconscan/runs/repro"

# ThermoRawFileParser (only needed to convert DDA .raw -> .mzML; ISD reads .raw directly)
TRFP="C:/Program Files/OpenMS-3.0.0-pre-HEAD-2023-06-17/share/OpenMS/THIRDPARTY/ThermoRawFileParser/ThermoRawFileParser.exe"

# ---- helpers -------------------------------------------------------------
# proteoform IDs at 1% FDR (and PTM-bearing = Full Sequence carries a modification)
count_ids(){ local psm; psm=$(find "$1" -name AllProteoforms.psmtsv 2>/dev/null | head -1)
  [ -z "$psm" ] && { echo "ERR ERR"; return; }
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="QValue")q=i;if($i=="Decoy/Contaminant/Target")t=i;if($i=="Full Sequence")fs=i}}
              NR>1 && $q+0<=0.01 && ($t=="T"||$t=="N"){tot++; if(index($fs,"[")>0)ptm++}
              END{print (tot+0)" "(ptm+0)}' "$psm"; }
