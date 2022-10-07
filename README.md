# GWAS and METAGWAS
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
## Meta-analysis
To increase the power of detection of significative signals, it is desirable to analyze multiple cohorts at once by running a Meta-GWAS. There are several tools to do it, in my case, I chose METAL.
1. **SCHEME** - To specify the type of analysis. *SampleSize* is the default option, using the effect, p-value and N of the study. It is recommendable when the analysis of the cohorts has been done with different tools. *STDERR* is the classical approach, using the effect estimations and the standard error.
2. **TRACKPOSITIONS** -  METAL checks if chromosome and position of a variant match across studies.
3. **GENOMICCONTROL** - To correct test statistics to account for small amounts of population stratification or unaccounted for relatedness
```
# SET PARAMETERS TO ANALYZE
SCHEME SAMPLESIZE/STDERR
TRACKPOSITIONS  ON
GENOMICCONTROL ON

#SET THE NAMES OF HEADERS THAT WE HAVE IN OUR FILES
MARKER -
POS -
WEIGHT -
ALLELE - -
FREQ -
EFFECT -
STDERR -
PVALUE -

# LOAD THE FILES
PROCESS $FILES

# CARRY OUT THE ANALYSIS
OUTFILE $OUTPUT .tbl
ANALYZE
ANALYZE HETEROGENEITY
```
![Manhattan Plot](https://upload.wikimedia.org/wikipedia/commons/4/4f/Manhattan_plot_from_a_GWAS_of_kidney_stone_disease.png)
