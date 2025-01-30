#!/bin/bash
echo "Reconstructing large files..."
cat app/out/0129b_part_* > app/out/0129b.csv
echo "Reconstruction complete."

