BASE=BASENAME

GENDER=temp/update_sex.txt
HIGH_LD_b37=help_files/high_LD_regions_b37.txt

SUBJ_FMISS_TH=0.03
FIRST_SNP_FMISS_TH=0.1
SNP_FMISS_TH=0.02
MAF_TH=0.005
HWE_TH=1e-06
HET_SD_TH=4
DUPL_TH=0.9

TEMP=./temp
HELP=./helper_scripts
RESULT=./result
SUMMARY=${RESULT}/${BASE}_QC_summary.txt

################################################################################
##                      Start writing to summary file                         ##
################################################################################
echo -e "### QC SUMMARY FOR PRESENT ###\n" > ${SUMMARY}
#echo "Input file: ${INPUT_PED_MAP}" >> ${SUMMARY}
echo "Affymetrix chip: Axiom_GW_Hu_SNP" >> ${SUMMARY}
echo "These files are on build 37 and + strand" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}.fam | awk '{print $1}'`
echo "Num subjects: ${NSAMPL}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}.bim | awk '{print $2}'`
echo -e "Num SNPs: ${NSNP}\n" >> ${SUMMARY}

################################################################################
##                      Check all sample IDs unique                           ##
################################################################################
echo "The following two output should be the same"
wc -l ${TEMP}/${BASE}.fam
cut -d' ' -f2 ${TEMP}/${BASE}.fam | sort -u | wc -l
echo "And this output should be empty"
cut -d' ' -f2 ${TEMP}/${BASE}.fam | sort | uniq -d
echo "---"

################################################################################
##                              Update FID                                    ##
################################################################################
module load cesga/2020 gcc/system plink/2.00a2.3
awk '{print $1, $2, $2, $2}' ${TEMP}/${BASE}.fam > ${TEMP}/${BASE}_update_ids.txt
plink2 --bfile ${TEMP}/${BASE} --update-ids ${TEMP}/${BASE}_update_ids.txt --make-bed --out ${TEMP}/${BASE}_raw_new_id

################################################################################
##                              Update sex                                    ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_raw_new_id --update-sex ${GENDER} --make-bed --out ${TEMP}/${BASE}_sex

################################################################################
##                       First exclusion of SNP missingness                   ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_sex --missing --out ${TEMP}/${BASE}_first_snp_miss
Rscript ${HELP}/lmiss_fail.R ${FIRST_SNP_FMISS_TH} ${TEMP}/${BASE}_first_snp_miss.vmiss ${TEMP}/${BASE}_first_fail_lmiss.txt
plink2 --bfile ${TEMP}/${BASE}_sex --exclude ${TEMP}/${BASE}_first_fail_lmiss.txt --make-bed --out ${TEMP}/${BASE}_first_snps_checked

## -- Write summary -- ##
echo "### First SNP call rate ###" >> ${SUMMARY}
echo "Threshold used: ${FIRST_SNP_FMISS_TH}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_first_fail_lmiss.txt | awk '{print $1}'`
echo -e "Num of SNPs failing: ${NSNP}\n" >> ${SUMMARY}

## Copy fail file
cp ${TEMP}/${BASE}_first_fail_lmiss.txt ${RESULT}

################################################################################
#                           Subject call rate                                 ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_first_snps_checked --missing --out ${TEMP}/${BASE}_subj_miss
tail -n +2 ${TEMP}/${BASE}_subj_miss.smiss | awk -v th=${SUBJ_FMISS_TH} '$6 > th {print $1, $2}' > ${TEMP}/${BASE}_fail_subj_callrate.txt

## -- Write summary -- ##
echo "### Sample call rate ###" >> ${SUMMARY}
echo "Threshold used: ${SUBJ_FMISS_TH}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_fail_subj_callrate.txt | awk '{print $1}'`
echo -e "Num of samples failing call rate: ${NSAMPL}\n" >> ${SUMMARY}

## Copy fail file
cp ${TEMP}/${BASE}_fail_subj_callrate.txt ${RESULT}

################################################################################
#                               Check sex                                     ##
################################################################################
## Independent markers
plink2 --bfile ${TEMP}/${BASE}_first_snps_checked --maf 0.1 --exclude range ${HIGH_LD_b37} --indep-pairwise 50 5 0.2 --out ${TEMP}/indepSNP_sexcheck

module load plink/1.9b5
plink --bfile ${TEMP}/${BASE}_first_snps_checked --extract ${TEMP}/indepSNP_sexcheck.prune.in --check-sex --out ${TEMP}/${BASE}_sexcheck1
plink --bfile ${TEMP}/${BASE}_first_snps_checked --extract ${TEMP}/indepSNP_sexcheck.prune.in --check-sex 0.3 0.8 --out ${TEMP}/${BASE}_sexcheck2
Rscript ${HELP}/sexcheck_fail.R ${TEMP}/${BASE}_sexcheck2.sexcheck ${TEMP}/${BASE}_fail_sex_shifted.txt ${TEMP}/${BASE}_fail_sex_uncertain.txt

# Exclude all
cat ${TEMP}/${BASE}_fail_sex_shifted.txt ${TEMP}/${BASE}_fail_sex_uncertain.txt > ${TEMP}/${BASE}_fail_sex.txt

# Write summary
echo "### Check sex ###" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_fail_sex.txt | awk '{print $1}'`
echo "Num of samples failing sex check: ${NSAMPL}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_fail_sex_shifted.txt | awk '{print $1}'`
echo "Num of samples with sex shifted: ${NSAMPL}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_fail_sex_uncertain.txt | awk '{print $1}'`
echo -e "Num of samples with uncertain sex: ${NSAMPL}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_sex_shifted.txt ${RESULT}
cp ${TEMP}/${BASE}_fail_sex_uncertain.txt ${RESULT}
cp ${TEMP}/${BASE}_fail_sex.txt ${RESULT}

################################################################################
##                       First exclusion of subjects                          ##
################################################################################
## without filtering out samples with wrong sex
cat ${TEMP}/${BASE}_fail_subj_callrate.txt ${TEMP}/${BASE}_fail_sex.txt | sort -u > ${TEMP}/${BASE}_first_subj_exclusion.txt
## filtering out those samples with wrong sex
cat ${TEMP}/${BASE}_fail_subj_callrate.txt ${TEMP}/${BASE}_fail_sex.txt | sort -u > ${TEMP}/${BASE}_first_subj_exclusion.txt
module load cesga/2020 gcc/system plink/2.00a2.3
plink2 --bfile ${TEMP}/${BASE}_first_snps_checked --remove ${TEMP}/${BASE}_first_subj_exclusion.txt --make-bed --out ${TEMP}/${BASE}_firstsubj_checked

## -- Write summary -- ##
echo "### First exclusion of subjects ###" >> ${SUMMARY}
echo "Excluding subjects failing call rate or sex check" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_first_subj_exclusion.txt | awk '{print $1}'`
echo "Num of samples failing: ${NSAMPL}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_firstsubj_checked.fam | awk '{print $1}'`
echo "Remaining subjects: ${NSAMPL}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_firstsubj_checked.bim | awk '{print $1}'`
echo -e "Remaining SNPs: ${NSNP}\n" >> ${SUMMARY}

## Copy fail file
cp ${TEMP}/${BASE}_first_subj_exclusion.txt ${RESULT}

################################################################################
##      Exclude true duplicated SNPS. Exclude the one with worst FMISS        ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_firstsubj_checked --missing --out ${TEMP}/${BASE}_snp_miss_dupl
Rscript ${HELP}/exclude_duplicate_snps.R ${TEMP}/${BASE}_firstsubj_checked.bim ${TEMP}/${BASE}_snp_miss_dupl.vmiss ${TEMP}/${BASE}_fail_dupl_snps.txt
plink2 --bfile ${TEMP}/${BASE}_firstsubj_checked --exclude ${TEMP}/${BASE}_fail_dupl_snps.txt --make-bed --out ${TEMP}/${BASE}_no_dupl_snps

# Write summary
echo "### Exclude duplicated SNPs ###" >> ${SUMMARY}
echo "When markers have the same chromosome position and the same alleles, exclude the ones with highest missingness" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_fail_dupl_snps.txt | awk '{print $1}'`
echo "Num of SNPs excluded: ${NSNP}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_no_dupl_snps.fam | awk '{print $1}'`
echo "Remaining subjects: ${NSAMPL}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_no_dupl_snps.bim | awk '{print $1}'`
echo -e "Remaining SNPs: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_dupl_snps.txt ${RESULT}

################################################################################
## 				mac=0 / monomorphic			      ##
################################################################################

plink2 --bfile ${TEMP}/${BASE}_no_dupl_snps --freq --out ${TEMP}/${BASE}_mac
Rscript ${HELP}/mac_fail.R ${TEMP}/${BASE}_mac.afreq ${TEMP}/${BASE}_fail_mac.txt
plink2 --bfile ${TEMP}/${BASE}_no_dupl_snps --exclude ${TEMP}/${BASE}_fail_mac.txt --make-bed --out ${TEMP}/${BASE}_mac_checked

# Write summary
echo "### Exclude monomorphic SNPs ###" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_fail_mac.txt | awk '{print $1}'`
echo "Num of SNPs excluded: ${NSNP}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_mac_checked.fam | awk '{print $1}'`
echo "Remaining subjects: ${NSAMPL}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_mac_checked.bim | awk '{print $1}'`
echo -e "Remaining SNPs: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_mac.txt ${RESULT}

################################################################################
# 				Minor allele freq			      ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_mac_checked --freq --out ${TEMP}/${BASE}_frq
Rscript ${HELP}/maf_fail.R ${MAF_TH} ${TEMP}/${BASE}_frq.afreq ${TEMP}/${BASE}_fail_maf.txt

# Write summary
echo "### Minor allele frequency ###" >> ${SUMMARY}
echo "Threshold used: ${MAF_TH}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_fail_maf.txt | awk '{print $1}'`
echo -e "Num of SNPs failing: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_maf.txt ${RESULT}

################################################################################
## 				SNP missingness				      ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_mac_checked --missing --out ${TEMP}/${BASE}_snp_miss
Rscript ${HELP}/lmiss_fail.R ${SNP_FMISS_TH} ${TEMP}/${BASE}_snp_miss.vmiss ${TEMP}/${BASE}_fail_lmiss.txt

# Write summary
echo "### SNP call rate ###" >> ${SUMMARY}
echo "Threshold used: ${SNP_FMISS_TH}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_fail_lmiss.txt | awk '{print $1}'`
echo -e "Num of SNPs failing: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_lmiss.txt ${RESULT}

################################################################################
## 					HWE				      ##
################################################################################
plink2 --bfile ${TEMP}/${BASE}_mac_checked --hardy --out ${TEMP}/${BASE}_hwe
Rscript ${HELP}/hwe_fail.R ${HWE_TH} ${TEMP}/${BASE}_hwe.hardy ${TEMP}/${BASE}_fail_hwe.txt

# Write summary
echo "### HWE ###" >> ${SUMMARY}
echo "Threshold used: ${HWE_TH}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_fail_hwe.txt | awk '{print $1}'`
echo -e "Num of SNPs failing: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_fail_hwe.txt ${RESULT}

################################################################################
## 				Exclusion of SNPs			      ##
################################################################################
cat ${TEMP}/${BASE}_fail_maf.txt ${TEMP}/${BASE}_fail_lmiss.txt ${TEMP}/${BASE}_fail_hwe.txt | sort -u > ${TEMP}/${BASE}_snp_exclusion.txt
plink2 --bfile ${TEMP}/${BASE}_mac_checked --exclude ${TEMP}/${BASE}_snp_exclusion.txt --make-bed --out ${TEMP}/${BASE}_snps_checked

# Write summary
echo "### Exclusion of SNPs ###" >> ${SUMMARY}
echo "Excluding SNPs failing MAF or call rate or HWE" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_snp_exclusion.txt | awk '{print $1}'`
echo "Num of SNPs failing: ${NSNP}" >> ${SUMMARY}
NSAMPL=`wc -l ${TEMP}/${BASE}_snps_checked.fam | awk '{print $1}'`
echo "Remaining subjects: ${NSAMPL}" >> ${SUMMARY}
NSNP=`wc -l ${TEMP}/${BASE}_snps_checked.bim | awk '{print $1}'`
echo -e "Remaining SNPs: ${NSNP}\n" >> ${SUMMARY}

# Copy fail file
cp ${TEMP}/${BASE}_snp_exclusion.txt ${RESULT}

################################################################################
## 			Autosomal indep SNPs				      ##
################################################################################
# Exclude non-autosomal
awk '($1 < 1) || ($1 > 22) {print $2}' ${TEMP}/${BASE}_snps_checked.bim > ${TEMP}/autosomeexcludes.txt
plink2 --bfile ${TEMP}/${BASE}_snps_checked --exclude ${TEMP}/autosomeexcludes.txt --make-bed --out ${TEMP}/${BASE}_autosomal

# Indep markers
plink2 --bfile ${TEMP}/${BASE}_autosomal --maf 0.1 --exclude range ${HIGH_LD_b37} --indep-pairwise 50 5 0.2 --out ${TEMP}/${BASE}_autosomal_indepSNP

################################################################################
# HapMap PCA for heterozygosity check
################################################################################
# Update to b37 forward strand using Will's script and files
# Copy strand b37 file
cp ./help_files/${STRAND_FILE} ${TEMP}

# Run Will's script
${HELP}/update_build_new.sh ${TEMP}/${BASE}_snps_checked ${TEMP}/${STRAND_FILE} ${TEMP}/${BASE}_hetcheck_b37fwd


# select only individuals from CEU, JPT, CHB, YRI
fam <- read.table('./help_files/HapMapIII_CGRCh37.fam')
pop <- read.delim("../reference/relationships_w_pops_121708.txt")
mm <- match(fam[,2],pop[,2])
fam[,6] <- as.character(pop$population[mm])
fam2 <- fam[fam[,6] %in% c('CEU')]


# create plink file with reduced populations
awk '{print $1, $2}' help_files/fam_CEU_JPT_CHB_YRI.fam > help_files/inds_CEU_JPT_CHB_YRI.txt

plink2 --bfile ./help_files/HapMapIII_CGRCh37 --keep help_files/inds_CEU_JPT_CHB_YRI.txt --make-bed --out ./help_files/HapMapIII_CGRCh37_reduced
