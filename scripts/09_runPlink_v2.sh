# run plink2 for the GWAS ------------------------s
# tool in CIIM 
/vol/projects/CIIM/resources/tools/plink2

# H1N1 
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_pheno_H1N1.txt \
  --glm allow-no-covars \
  --1

# H3N2
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_pheno_H3N2.txt \
  --glm allow-no-covars \
  --1

# B
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_pheno_B.txt \
  --glm allow-no-covars \
  --1

# extract sig. SNPs --------------------------
# because the vcf file have some preflix in the SNPs, so could not use the bcftools, we use the defaut linux command
grep -w -e CHROM -f /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/qtl_sigSNPs_H3N2_onlyRsid.txt iMED_vcf.vcf >| /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/processedDat/qtl_sigSNPs_H3N2_genotype.txt 

# run plink2 - QTL with abTiter at T4 ------------------------s
# H1N1 
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_H1N1_T4log2.txt \
  --glm allow-no-covars \
  --1

# H3N2
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_H3N2_T4log2.txt \
  --glm allow-no-covars \
  --1
  
# B 
/vol/projects/CIIM/resources/tools/plink2 \
  --vcf /vol/projects/CIIM/Influenza/iMED/genotype/genotypes_combined_new/qtl_mapping/iMED_vcf.vcf \
  --pheno iid-only /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/sendOut/iMED_B_T4log2.txt \
  --glm allow-no-covars \
  --1