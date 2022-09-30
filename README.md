# GWAS
Scripts to run Genome-Wide Association Studies and COX-GWAS.

## 

```
plink2 --bfile $GENODATA --maf 0.05 --geno 0.02 --mind 0.02 --hwe 1e-6 --make-bed --out $OUTPUT_FILTERED
```

```
plink2 --bfile $OUTPUT_FILTERED --pheno $PHENO --keep $IDS --covar $COVS --glm hide-covar --freq --threads 64 --out $OUTPUT
```
