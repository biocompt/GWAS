MAIN="path/to/main_folder"
RAW="path/to/rawdata"
UNFILTERED="path/to/unfiltered"
FILTERED="path/to/filtered"
CONCATENATED="path/to/concatenated"
SCRIPTS="path/to/scripts_folder"

## 1) Extract each sample from the zip file downloaded
## We will create separately scripts for each chromosome, so we can run all of them in one round
ls $RAW > $SCRIPTS/raw_chr.txt
mkdir $SCRIPTS/UNZIP
for line in $SCRIPTS/raw_chr.txt
do
echo "!/bin/bash" > $SCRIPTS/UNZIP/${line::-4}.sh
echo "unzip -P 'PASSWORD' $line -d ${UNFILTERED}" >> $SCRIPTS/UNZIP/${line::-4}.sh
done

for file in $SCRIPTS/UNZIP/*.sh
do
sbatch $file
done

echo "You have unzipped all the files"

## 2) Change name of the file to describe it
cd $UNFILTERED
# Move info files to a different folder
mkdir ../info
mv *.info ../info

for file in *.vcf.gz
do
mv $file $(basename $file .dose.vcf.gz)_[DESCRIPTION_OF_THE_FILE].dose.vcf.gz
done

## 3) Indexing tab-delimited genome position
ls $UNFILTERED > $SCRIPTS/unzipped_chr.txt
mkdir $SCRIPTS/TABIX
for line in $SCRIPTS/unzipped_chr.txt
do
echo "!/bin/bash" > $SCRIPTS/TABIX/${line::-7}.sh
echo "tabix -p vcf $file" >> $SCRIPTS/TABIX/${line::-7}.sh
done

for file in $SCRIPTS/TABIX/*.sh
do
sbatch $file
done

## 4) Filter by INFO imputation score
mkdir $SCRIPTS/FILTERED 
for line in $SCRIPTS/unzipped_chr.txt
do
echo "!/bin/bash" > $SCRIPTS/TABIX/${line::-7}.sh
echo "bcftools view -i 'MIN(INFO/R2)>0.3' $line --threads 8 | bgzip -c > $FILTERED/${line::-7}.INFO_0.3.vcf.gz" >> $SCRIPTS/FILTERED/${line::-7}.sh
done

for file in $SCRIPTS/FILTERED/*.sh
do
sbatch $file
done

## 5) Indexing filtered chromosomes
## Repeat step 3

## 6) Concatenated all chromosomes in one file
bcftools concat chr*.vcf.gz -o $CONCATENATED/OUTPUT
