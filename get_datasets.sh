#!/usr/bin/env bash
set -e

DST_ROOT=${1:-$HOME/dihiggs/datasets}
HB_URL="https://gitlab.com/higgsbounds/hbdataset.git"
HS_URL="https://gitlab.com/higgsbounds/hsdataset.git"

mkdir -p "$DST_ROOT"

if [ ! -d "$DST_ROOT/HBDataset" ]; then
  echo "Cloning HBDataset..."
  git clone --depth 1 "$HB_URL" "$DST_ROOT/HBDataset"
else
  echo "HBDataset already present, pulling updates..."
  cd "$DST_ROOT/HBDataset"
  git pull
  cd -
fi

if [ ! -d "$DST_ROOT/HSDataset" ]; then
  echo "Cloning HSDataset..."
  git clone --depth 1 "$HS_URL" "$DST_ROOT/HSDataset"
else
  echo "HSDataset already present, pulling updates..."
  cd "$DST_ROOT/HSDataset"
  git pull
  cd -
fi

echo "Datasets ready in $DST_ROOT"
