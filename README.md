# GWAS
Scripts to run Genome-Wide Association Studies and COX-GWAS.

## Filtering the genotype
It is important to filter the initial file to ensure a reliable analysis.
1. **--maf**: 0.05 -> Filters out all variants with minor allele frequency below 5%. 
2. **--geno**: 0.02 -> Filters out all variants with missing call rate higher than 2%.
3. **--mind**: 0.02 -> Filters out all samples with missing call rate higher than 2%.
4. **--hwe**: 1e-6 -> Filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below 1e-6
```
plink2 --bfile $GENODATA --maf 0.05 --geno 0.02 --mind 0.02 --hwe 1e-6 \
  --make-bed --out $OUTPUT_FILTERED
```
## Run the analysis
We set the desired parameters to our study. As an option, we can give a pheno file with multiple phenos, one for colum. Usually, we are just interest in obtain the "ADD" values in the output, so we can run fast the analysis and obtaining a smaller output by setting the option hide-covar.
```
plink2 --bfile $OUTPUT_FILTERED --pheno $PHENO --keep $IDS --covar $COVS \
  --glm hide-covar cols=chrom,pos,ref,alt,ax,a1freq,nobs,beta,se,p --out $OUTPUT
```
