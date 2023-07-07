#!/bin/bash

inputfile=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/input/iMED_updatedID
outfile=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/output
models=/vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/omicsPred_models

for datType in Metabolon Nightingale Olink RNAseq SomaScan; 
do 
  mkdir $outfile/$datType
  for i in $(ls $models/$datType); do
  /vol/projects/CIIM/resources/tools/plink2 --bad-freqs \
    --bfile $inputfile \
    --score $models/$datType/$i 1 4 6 header list-variants cols=scoresums \
    --out $outfile/$datType/$i ; 
  done
done
