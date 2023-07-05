#!/bin/bash
cp /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/input

cd /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/omicsPred/input

# Step 1: bcftools query -f "%ID" iMED_vcf.vcf # List of SNP IDs per SNP. Save in a file
/vol/projects/CIIM/resources/tools/bcftools query -f '%ID\n' iMED_vcf.vcf > iMED_IDs.txt

# Step 2: Take file from step 1, slpit on % (to get only the rsid without Chro position). 'Paste' together with file from step 1
awk '{print}' iMED_IDs.txt | sed s/.*%// > iMED_rsid.txt
paste iMED_IDs.txt iMED_rsid.txt > iMED_ID_rsid.txt

# Step 3: generate .fam, .bim. bed files from .vcf file with new rsID (using .txt file from step 2) and without rewrite the .vcf files
/vol/projects/CIIM/resources/tools/plink2 --vcf iMED_vcf.vcf \
  --update-name iMED_ID_rsid.txt --vcf-idspace-to _ \
  --make-bed --out iMED_updatedID
