#!/bin/bash

inputfile=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/input/iMED_updatedID
outfile=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/Olink_output/
models=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/Olink_OMICSPRED

for i in $(ls $models); do
/vol/projects/CIIM/resources/tools/plink2 --bad-freqs \
  --bfile $inputfile \
  --score $models/$i 1 4 6 header list-variants cols=scoresums \
  --out $outfile$i ; 
done