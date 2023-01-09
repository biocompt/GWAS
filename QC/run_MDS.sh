##############################################
###########    GENOTYPING DATA    ############
########### pre-imputation MDS QC ############
###########                	  ############
##############################################

KG="1kG/1kG"
GENO="geno/REQUITE_TOPMed_2022_GENO_0.02_MIND_0.02_MAF_0.05_HWE_HET_FOUNDERS"
PCA="pca"

module load cesga/2020 gcc/system plink/2.00a2.3
## QC on 1000 Genomes data.
# Remove variants based on missing genotype data.
plink2 --bfile $KG --geno 0.2 --allow-no-sex --make-bed --out ${KG}_MDS

# Remove individuals based on missing genotype data.
plink2 --bfile ${KG}_MDS --mind 0.2 --allow-no-sex --make-bed --out ${KG}_MDS2

# Remove variants based on missing genotype data.
plink2 --bfile ${KG}_MDS2 --geno 0.02 --allow-no-sex --make-bed --out ${KG}_MDS3

# Remove individuals based on missing genotype data.
plink2 --bfile ${KG}_MDS3 --mind 0.02 --allow-no-sex --make-bed --out ${KG}_MDS4

# Remove variants based on MAF.
plink2 --bfile ${KG}_MDS4 --maf 0.05 --allow-no-sex --make-bed --out ${KG}_MDS5

module load plink/1.9b5
#>>>>>>>>>>>>>>>>>>>>>>>>> DESDE AQUI SE PUEDE INICIAR CON EL 1kG_MDS5
#Como los ids de los snps estan en formato distinto en los 2 archivos
#cambiamos el de los 1kG al formato de nuestra cohorte
plink --bfile 1kG/1kG --update-name new_ids.txt --make-bed --out 1kG_updated

module load cesga/2020 gcc/system plink/2.00a2.3
# Extract the variants present in HapMap dataset from the 1000 genomes dataset.
awk '{print$2}' geno/4_COHORTES.bim > HapMap_SNPs.txt
plink2 --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6

module load plink/1.9b5
# Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
#----> ELIMINO ambos SNPs DUPLICADOS para evitar errores | 4:52938243, 9:135770300, 9:135770347, 9:136135237
plink --bfile geno/4_COHORTES --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out HapMap_MDS
# The datasets now contain the exact same variants.

module load cesga/2020 gcc/system plink/2.00a2.3
## The datasets must have the same build. Change the build 1000 Genomes data build.
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
# buildhapmap.txt contains one SNP-id and physical position per line.
plink2 --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7

# 1kG_MDS7 and HapMap_MDS now have the same build.
###################################################################

## Merge the HapMap and 1000 Genomes data sets
# Prior to merging 1000 Genomes data with the HapMap data we want to make sure that the files are mergeable, for this we conduct 3 steps:
# 1) Make sure the reference genome is similar in the HapMap and the 1000 Genomes Project datasets.
# 2) Resolve strand issues.
# 3) Remove the SNPs which after the previous two steps still differ between datasets.

# The following steps are maybe quite technical in terms of commands, but we just compare the two data sets and make sure they correspond.

module load plink/1.9b5
# 1) set reference genome
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
# The 1kG_MDS7 and the HapMap-adj have the same reference genome for all SNPs.
# This command will generate some warnings for impossible A1 allele assignment.
plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
sort 1kGMDS7_tmp HapMap-adj_tmp | uniq -u > all_differences.txt
# 2 differences between the files, some of these might be due to strand issues.

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
# Generates a file of 1 SNPs. These are the non-corresponding SNPs between the two files.
# Flip the 1 non-corresponding SNPs.
plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
# This file demonstrates that there are 84 differences between the files.

# 3) Remove problematic SNPs from HapMap and 1000 Genomes.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
# The command above generates a list of the 42 SNPs which caused the 84 differences between the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

module load cesga/2020 gcc/system plink/2.00a2.3
# Remove the 42 problematic SNPs from both datasets.
plink2 --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2
plink2 --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8

# Merge HapMap with 1000 Genomes Data.
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2

# Note, we are fully aware of the sample overlap between the HapMap and 1000 Genomes datasets. However, for the purpose of this tutorial this is not important.

## Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
# Using a set of pruned SNPs
plink --bfile HapMap_MDS2 --genome --out MDS_merge3
plink --bfile HapMap_MDS2 --read-genome MDS_merge3.genome --cluster --mds-plot 10 --out MDS_merge2

### MDS-plot

# Download the file with population information of the 1000 genomes dataset.
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# The file 20100804.ALL.panel contains population codes of the individuals of 1000 genomes.
# Extraemos el ID y la superpoblacion a la que pertenece
mv integrated_call_samples_v3.20130502.ALL.panel ALL_panel
cat ALL_panel | cut -f 1,3 > races.txt

# Create a racefile of your own data.
awk '{print$1,$2,"OWN"}' HapMap_MDS.fam > racefile_own.txt

# Concatenate racefiles.
cat races.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

# Generate population stratification plot.
Rscript MDS_merged.R
# The output file MDS.pdf demonstrates that our "owndata falls within the European group of the 1000 genomes data. Therefore, we do not have to remove subjects.
# For educational purposes however, we give scripts below to filter out population stratification outliers. Please execute the script below in order to generate the appropriate files for the next tutorial.

## Exclude ethnic outliers.
# Select individuals in HapMap data below cut-off thresholds. The cut-off levels are not fixed thresholds but have to be determined based on the visualization of the first two dimensions. To exclude ethnic outliers, the thresholds need to be set around the cluster of population of interest.
awk '{ if ($4 < -0.03 && $5 > 0.02) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2

module load cesga/2020 gcc/system plink/2.00a2.3
# Extract these individuals in HapMap data.
plink2 --bfile IBIS_v2_geno0.01_mind0.025_sexcheck_autosomes_maf0.05_hwe106_hh3sd_founders_pihat0.2 --make-bed --keep EUR_MDS_merge2 --out IBIS_v2_2021
#Start time: Wed Sep  1 22:30:04 2021
#32557 MiB RAM detected; reserving 16278 MiB for main workspace.
#Using up to 12 threads (change this with --threads).
#185 samples (29 females, 156 males; 185 founders) loaded from
#IBIS_v2_geno0.01_mind0.025_sexcheck_autosomes_maf0.05_hwe106_hh3sd_founders_pihat0.2.fam.
#351779 variants loaded from
#IBIS_v2_geno0.01_mind0.025_sexcheck_autosomes_maf0.05_hwe106_hh3sd_founders_pihat0.2.bim.
#Note: No phenotype data present.
#--keep: 182 samples remaining.
#182 samples (28 females, 154 males; 182 founders) remaining after main filters.

plink1.9 --bfile IBIS_v2_2021 --extract indepSNP.prune.in --genome --out IBIS_v2_2021
#351779 variants loaded from .bim file.
#182 people (154 males, 28 females) loaded from .fam.
#--extract: 123566 variants remaining.
#Using up to 11 threads (change this with --threads).
#Before main variant filters, 182 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.999035.
#123566 variants and 182 people pass filters and QC.
#Note: No phenotypes present.
#IBD calculations complete.
#Finished writing IBIS_v2_2021.genome .

plink1.9 --bfile IBIS_v2_2021 --read-genome IBIS_v2_2021.genome --cluster --mds-plot 10 --out IBIS_v2_2021_mds
#351779 variants loaded from .bim file.
#182 people (154 males, 28 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 182 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.999073.
#351779 variants and 182 people pass filters and QC.
#Note: No phenotypes present.
#Clustering... done.
#Cluster solution written to IBIS_v2_2021_mds.cluster1 ,
#IBIS_v2_2021_mds.cluster2 , and IBIS_v2_2021_mds.cluster3 .
#Performing multidimensional scaling analysis (SVD algorithm, 10 dimensions)... done.
#MDS solution written to IBIS_v2_2021_mds.mds .

# Change the format of the .mds file into a plink covariate file.
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' IBIS_v2_2021_mds.mds > covar_mds.txt
# The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.

##########################################################################################################################################################################

## CONGRATULATIONS you have succesfully controlled your data for population stratification!
