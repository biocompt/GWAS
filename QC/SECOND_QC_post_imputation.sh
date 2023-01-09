MAIN="path/to/main_folder"
RAW="path/to/rawdata"
UNFILTERED="path/to/unfiltered"
FILTERED="path/to/filtered"
CONCATENATED="path/to/concatenated"
SCRIPTS="path/to/scripts_folder"

## 1) Extract each sample from the zip file downloaded
cd $RAW
for i in {1..22}
do
unzip -P "password" chr_${i}.zip -d $UNFILTERED/ &
done

## 2) Change name of the file to describe it
cd $UNFILTERED
for file in *.vcf.gz
do
mv $file $(basename $file .dose.vcf.gz)_[DESCRIPTION_OF_THE_FILE].dose.vcf.gz
done

# Move info files to a different folder
mkdir ../info
mv *.info ../info

## 3) Indexing tab-delimited genome position
for file in *.vcf.gz
do
tabix -p vcf $file &
done

## 4) Filter by INFO imputation score
for file in *.vcf.gz
do
bcftools view -i 'MIN(INFO/R2)>0.3' chr${chr}_[DESCRIPTION_OF_THE_FILE].dose.vcf.gz --threads 8 | bgzip -c > $FILTERED/chr${chr}_[DESCRIPTION_OF_THE_FILE].dose.INFO0.3.vcf.gz
done

## 5) Indexing filtered chromosomes
cd $FILTERED
for file in *.vcf.gz
do
tabix -p vcf $file &
done

## 6) Concatenated all chromosomes in one file
bcftools concat chr*.vcf.gz -o $CONCATENATED/OUTPUT
