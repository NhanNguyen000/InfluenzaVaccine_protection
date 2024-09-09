# Generate SNPloc input file for magma annotate step

mkdir B

# Contains three cols: SNPid, chr, bp on hg19
cat B_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,1,2 | awk -v OFS='\t' '{print $3, $1, $2}'> B/snploc.txt

# Annotate SNPs
/vol/projects/CIIM/resources/tools/MAGMA/magma --annotate window=10,10 --snp-loc B/snploc.txt --gene-loc /vol/projects/CIIM/resources/tools/MAGMA/aux_files/NCBI38.gene.loc --out B/step1

# Gene level
echo -e 'SNP\tP' > B/pvals.txt
cat B_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,9 >> B/pvals.txt
/vol/projects/CIIM/resources/tools/MAGMA/magma --bfile /vol/projects/CIIM/resources/tools/MAGMA/aux_files/g1000_eur --gene-annot B/step1.genes.annot --pval B/pvals.txt N=160 --out B/genelevel

# Gene set level
/vol/projects/CIIM/resources/tools/MAGMA/magma --gene-results B/genelevel.genes.raw --set-annot pathways.txt --out B/genesetlevel




# H3N2
mkdir H3N2

# Contains three cols: SNPid, chr, bp on hg19
cat H3N2_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,1,2 | awk -v OFS='\t' '{print $3, $1, $2}'> H3N2/snploc.txt

# Annotate SNPs
/vol/projects/CIIM/resources/tools/MAGMA/magma --annotate window=10,10 --snp-loc H3N2/snploc.txt --gene-loc /vol/projects/CIIM/resources/tools/MAGMA/aux_files/NCBI38.gene.loc --out H3N2/step1

# Gene level
echo -e 'SNP\tP' > H3N2/pvals.txt
cat H3N2_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,9 >> H3N2/pvals.txt
/vol/projects/CIIM/resources/tools/MAGMA/magma --bfile /vol/projects/CIIM/resources/tools/MAGMA/aux_files/g1000_eur --gene-annot H3N2/step1.genes.annot --pval H3N2/pvals.txt N=160 --out H3N2/genelevel

# Gene set level
/vol/projects/CIIM/resources/tools/MAGMA/magma --gene-results H3N2/genelevel.genes.raw --set-annot pathways.txt --out H3N2/genesetlevel







# H1N1
mkdir H1N1

# Contains three cols: SNPid, chr, bp on hg19
cat H1N1_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,1,2 | awk -v OFS='\t' '{print $3, $1, $2}'> H1N1/snploc.txt

# Annotate SNPs
/vol/projects/CIIM/resources/tools/MAGMA/magma --annotate window=10,10 --snp-loc H1N1/snploc.txt --gene-loc /vol/projects/CIIM/resources/tools/MAGMA/aux_files/NCBI38.gene.loc --out H1N1/step1

# Gene level
echo -e 'SNP\tP' > H1N1/pvals.txt
cat H1N1_T4.txt | tr "%" "\t" | tr ":" "\t" | grep -v  'SNP' | cut -f5,9 >> H1N1/pvals.txt
/vol/projects/CIIM/resources/tools/MAGMA/magma --bfile /vol/projects/CIIM/resources/tools/MAGMA/aux_files/g1000_eur --gene-annot H1N1/step1.genes.annot --pval H1N1/pvals.txt N=160 --out H1N1/genelevel

# Gene set level
/vol/projects/CIIM/resources/tools/MAGMA/magma --gene-results H1N1/genelevel.genes.raw --set-annot pathways.txt --out H1N1/genesetlevel
