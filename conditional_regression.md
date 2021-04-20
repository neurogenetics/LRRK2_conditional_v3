## Logistic regression model with LRRK2 coding variants as covariates

```
# Julie
# April 2021
```

### IPDGC
```
cd /PATH/TO/Julie/Julie_LRRK2_Condi
mkdir LRRK2_var_regression
cd LRRK2_var_regression

# Recode the genotypes of interest as single allele dosage numbers 
# Will use some of these as covariates in the GWAS
module load plink
plink --bfile /PATH/TO/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--snps 12:40614434,12:40629436,12:40657700,12:40671989,12:40702911,12:40707778,12:40707861,12:40713899,12:40734202,12:40740686 \
--recodeA --out LRRK2_rare_variant_selection

module load R
R
require(dplyr)
require(data.table)
vars <- fread("LRRK2_rare_variant_selection.raw",header=T)
vars$IID <- vars$PAT <- vars$MAT <- vars$SEX <- vars$PHENOTYPE <- NULL

# Add the variant info to the NORMAL covariate files, to be used as covariates
file.names <- dir("/PATH/TO/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES", pattern="LRRK2_condi_covariates_NORMAL.*.txt") 
for(file in file.names) {  
  cov <- fread(paste("/PATH/TO/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/", file,sep=""))
  Mrg <- merge(cov,vars,by='FID')
  write.table(Mrg, file=paste("vars",file,sep="_"), quote=FALSE,row.names=F,sep="\t")
}
q()
n

# NORMAL
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5 NORMAL_LRRK2_GWAS
# G2019S and rs76904798
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40614434_T,12:40734202_A GSRS_LRRK2_GWAS
# G2019S and N2081D
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40734202_A,12:40740686_G GSND_LRRK2_GWAS
# G2019S only
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40734202_A GS_LRRK2_GWAS
# rs76904798 only
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40614434_T RS_LRRK2_GWAS
# N2081D only
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40740686_G ND_LRRK2_GWAS
# LRRK2 PD-linked
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40657700_G,12:40702911_A,12:40713899_C,12:40734202_A,12:40740686_G PDlinked_LRRK2_GWAS
# LRRK2 rare
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=1:00:00 IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5,12:40629436_C,12:40657700_G,12:40671989_G,12:40702911_A,12:40707778_A,12:40707861_T,12:40713899_C,12:40734202_A,12:40740686_G Rare_LRRK2_GWAS
```
```
# This is IPDGC_LRRK2_GWAS.sh

#!/bin/bash
# sh IPDGC_LRRK2_GWAS.sh AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5 NORMAL_LRRK2_GWAS
COVARIATES=$1
OUTFILE=$2

cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/
mkdir ${OUTFILE}
cd ${OUTFILE}

module load plink/2.0-dev-20191128
cat /PATH/TO/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
	plink2 --bfile /PATH/TO/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--snps 12:40614434,12:40629436,12:40657700,12:40671989,12:40702911,12:40707778,12:40707861,12:40713899,12:40713901,12:40734202,12:40740686,12:40758652,12:123326598 \
	--keep /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.$line.txt \
	--pheno-name PHENO_PLINK --covar-variance-standardize \
	--pheno /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar-name $COVARIATES \
	--out $OUTFILE.$line
done

# EXCEPTIONS => VANCE + MF no age...
plink2 --bfile /PATH/TO/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snps 12:40614434,12:40629436,12:40657700,12:40671989,12:40702911,12:40707778,12:40707861,12:40713899,12:40713901,12:40734202,12:40740686,12:40758652,12:123326598 \
--keep /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.VANCE.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar-name $COVARIATES \
--out $OUTFILE.VANCE

plink2 --bfile /PATH/TO/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snps 12:40614434,12:40629436,12:40657700,12:40671989,12:40702911,12:40707778,12:40707861,12:40713899,12:40713901,12:40734202,12:40740686,12:40758652,12:123326598 \
--keep /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.MF.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.MF.txt \
--covar /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/vars_LRRK2_condi_covariates_NORMAL.MF.txt \
--covar-name $COVARIATES \
--out $OUTFILE.MF
```
```
# Reformat the GWAS results
sh reformat_and_organize.sh NORMAL_LRRK2_GWAS
sh reformat_and_organize.sh GSRS_LRRK2_GWAS
sh reformat_and_organize.sh GSND_LRRK2_GWAS
sh reformat_and_organize.sh GS_LRRK2_GWAS
sh reformat_and_organize.sh RS_LRRK2_GWAS
sh reformat_and_organize.sh ND_LRRK2_GWAS
sh reformat_and_organize.sh PDlinked_LRRK2_GWAS
sh reformat_and_organize.sh Rare_LRRK2_GWAS

# This is reformat_and_organize.sh

#!/bin/bash
# sh reformat_and_organize.sh NORMAL_LRRK2_GWAS
GWAS_TYPE=$1

# Loop over to reformat the .hybrid files into .txt files
cd /PATH/TOJulie/Julie_LRRK2_Condi/LRRK2_var_regression/${GWAS_TYPE}
cat /PATH/TO/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
  Rscript --vanilla /PATH/TO/Julie/Julie_LRRK2_Condi/reformat_IPDGC.R $GWAS_TYPE.$line.PHENO_PLINK.glm.logistic.hybrid
done

# Organize the files 
mkdir prep_files
mv *.log prep_files
mv *.hybrid prep_files
```

### UK Biobank
```
## Determine LRRK2 carrier status
cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression
mkdir UKB
cd UKB

module load plink/2.3-alpha
plink2 --bgen /PATH/TO/UKBIOBANK/IMPUTED_DATA/ukb_imp_chr12_v3.bgen \
--snps rs76904798,rs33995463, rs7308720, rs10878307, rs7133914, rs35507033, rs33958906, rs35303786, rs34637584, rs33995883 \
--make-bed --sample /PATH/TO/UKBIOBANK/IMPUTED_DATA/ukb33601_imp_chr1_v3_s487395.sample --out LRRK2_rare_variant_selection_UKB

module load plink
plink --bfile LRRK2_rare_variant_selection_UKB --recodeA --out LRRK2_rare_variant_selection_UKB2

# Add the variant info to the NORMAL cov files
cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB
module load R
R
require(dplyr)
require(data.table)
PD <- fread("/PATH/TO/Julie/Julie_LRRK2_Condi/UKB_GWAS/COV_UKB_PD_cases_control_over60.txt",header=T)
Proxy <- fread("/PATH/TO/Julie/Julie_LRRK2_Condi/UKB_GWAS/COV_UKB_Proxy_cases_control_over60.txt",header=T)
vars <- fread("LRRK2_rare_variant_selection_UKB2.raw",header=T)
vars$IID <- vars$PAT <- vars$MAT <- vars$SEX <- vars$PHENOTYPE <- NULL

# These are repeats in the cov file already
vars$rs76904798_T <- vars$rs34637584_A <- vars$rs33995883_G <- NULL

Mrg_PD <- merge(PD,vars,by='FID')
write.table(Mrg_PD, file="vars_COV_UKB_PD_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
Mrg_Proxy <- merge(Proxy,vars,by='FID')
write.table(Mrg_Proxy, file="vars_COV_UKB_Proxy_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")

# Perform GWAS
module load plink/2.0-dev-20191128

# NORMAL
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 NORMAL_LRRK2_GWAS
# G2019S and rs76904798
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs76904798_T,rs34637584_A GSRS_LRRK2_GWAS
# G2019S and N2081D
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs34637584_A,rs33995883_G GSND_LRRK2_GWAS
# G2019S only
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs34637584_A GS_LRRK2_GWAS
# rs76904798 only
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs76904798_T RS_LRRK2_GWAS
# N2081D only
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs33995883_G ND_LRRK2_GWAS
# LRRK2 PD-linked
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs7308720_G,rs7133914_A,rs35303786_C,rs34637584_A,rs33995883_G PDlinked_LRRK2_GWAS
# LRRK2 rare
sh UKB_LRRK2_GWAS.sh AGE_OF_RECRUIT,GENETIC_SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5,rs33995463_C,rs7308720_G,rs10878307_G,rs7133914_A,rs35507033_A,rs33958906_T,rs35303786_C,rs34637584_A,rs33995883_G Rare_LRRK2_GWAS
```
```
# This is UKB_LRRK2_GWAS.sh

#!/bin/bash
# sh UKB_LRRK2_GWAS.sh NORMAL_LRRK2_GWAS
COVARIATES=$1
OUTFILE=$2

cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB
plink2 --pfile /PATH/TO/Julie/Julie_LRRK2_Condi/UKB_GWAS/chr12.UKBB.EU.filtered_NEW \
--pheno-name STATUS --pheno vars_COV_UKB_PD_cases_control_over60.txt \
--covar vars_COV_UKB_PD_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snps rs76904798,rs33995463, rs7308720, rs10878307, rs7133914, rs35507033, rs33958906, rs35303786, rs11564148, rs34637584, rs33995883,rs3761863,rs10847864 \
--out $OUTFILE.UKBPD \
--covar-name $COVARIATES --covar-variance-standardize

plink2 --pfile /PATH/TO/Julie/Julie_LRRK2_Condi/UKB_GWAS/chr12.UKBB.EU.filtered_NEW \
--pheno-name STATUS --pheno vars_COV_UKB_Proxy_cases_control_over60.txt \
--covar vars_COV_UKB_Proxy_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snps rs76904798,rs33995463, rs7308720, rs10878307, rs7133914, rs35507033, rs33958906, rs35303786, rs11564148, rs34637584, rs33995883,rs3761863,rs10847864 \
--out $OUTFILE.UKBProxy \
--covar-name $COVARIATES --covar-variance-standardize
```
```
# Reformat plink2 GWAS output
cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB
ls *LRRK2_GWAS*.hybrid > GWAS_files.txt
module load python/3.6
cat GWAS_files.txt | while read line
do 
  # Filter the results by A1_FREQ (minor allele frequency) >=0.0001 --> to input.txt \
  awk '{ if($13 >= 0.0001) { print }}' $line > input.txt
  # Reformat the plink2 results for meta-analysis using python
  if [[ $line == *"PD"* ]]; then
    python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
    --outfile toMeta.${line%%.*}.UKBPD.txt --B-or-C B
  elif [[ $line == *"Proxy"* ]]; then
    python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
    --outfile toProxy.${line%%.*}.UKBProxy.txt --B-or-C B; fi
done

# Adjust Proxy cases to the same scale as PD cases
# Make the toProxy files into .csv files
module load R
R
require(data.table)
require(dplyr)
file.names <- dir("/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB", pattern="toProxy") 
for(file in file.names) {  
  data <- fread(file,header=T)
  newname <- file %>% gsub(pattern=".txt", replacement=".csv") %>% gsub(pattern="toProxy",replacement="toConvert")
  write.table(data, file=newname,quote=F,row.names=F,sep=",")
}
q()
n

# Convert proxies to "normal" files ready for meta-analysis
ls toConvert* > convert_files.txt
cat convert_files.txt | while read line
do
  pattern="toConvert"
  replacement="toMeta"
  newname="${line/$pattern/$replacement}"
  python /PATH/TO/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
  --infile $line --beta-proxy beta --se-proxy se --p-proxy P --outfile $newname
done

# Convert the "normalized" proxy .csv files back to .txt files
R
require(data.table)
require(dplyr)
file.names <- dir("/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB", pattern="toMeta.*UKBProxy.csv") 
for(file in file.names) {  
  data <- fread(file,header=T)
  newname <- file %>% gsub(pattern=".csv", replacement=".txt")
  write.table(data, file=newname,quote=F,row.names=F,sep="\t")
}
q()
n

# Reformat and move to the same directory as the IPDGC CHR12 files for meta-analysis
cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/UKB
ls toMeta*.txt > meta_files.txt
cat meta_files.txt | while read line
do
  # Get rid of the :N:N off of the markername so that it matches IPDGC
  sed -i -r 's/:[ACGT]+//g' $line
  pattern1="toMeta."
  pattern2=".UKB*txt"
  replacement=""
  newname1="${line/$pattern1/$replacement}"
  newname2="${newname1/$pattern2/$replacement}"
  if [[ $line == *"PD"* ]]; then
  cut -f1-7 $line | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/${newname2}/$newname2.UKBPD.txt
  elif [[ $line == *"Proxy"* ]]; then
  cut -f1,2,3,7,15,16,17 $line | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/${newname2}/$newname2.UKBProxy.txt; fi
done
```

### Meta-analysis in METAL
```
# Run meta-analysis in METAL
cd /PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression
mkdir METAL
cd METAL

sh /PATH/TO/Julie/METAL.sh "NORMAL" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/NORMAL_LRRK2_GWAS/NORMAL_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "GSRS" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/GSRS_LRRK2_GWAS/GSRS_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "GSND" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/GSND_LRRK2_GWAS/GSND_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "GS" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/GS_LRRK2_GWAS/GS_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "RS" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/RS_LRRK2_GWAS/RS_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "ND" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/ND_LRRK2_GWAS/ND_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "PDlinked" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/PDlinked_LRRK2_GWAS/PDlinked_LRRK2_GWAS.*.txt"
sh /PATH/TO/Julie/METAL.sh "Rare" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/Rare_LRRK2_GWAS/Rare_LRRK2_GWAS.*.txt"

# This is METAL.sh
#!/bin/bash
# sh METAL.sh "NORMAL" "/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/NORMAL_LRRK2_GWAS/NORMAL_LRRK2_GWAS.*.txt"
PREFIX=$1
PATTERN=$2

ls $PATTERN > ${PREFIX}_GWAS_FILES.txt
get_colname () {
colnumber=$1
file=$2
head -1 $file |cut -f$colnumber
}

make_metal_file () {
echo """
#../generic-metal/metal metalAll.txt
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile17.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON
"""

cat ${PREFIX}_GWAS_FILES.txt | while read line
do
	echo ""
	echo "# === DESCRIBE AND PROCESS INPUT FILE ==="
	echo MARKER $(get_colname 1 $line)
	echo ALLELE $(get_colname 3 $line) $(get_colname 2 $line)
	echo FREQ   $(get_colname 4 $line)
	echo EFFECT $(get_colname 5 $line)
	echo STDERR $(get_colname 6 $line)
	echo PVALUE $(get_colname 7 $line)
	echo PROCESS $line
done
echo """
OUTFILE ${PREFIX}_metal .tbl
ANALYZE HETEROGENEITY
QUIT
"""
}

module load metal
make_metal_file > ${PREFIX}_metal_template.txt
metal ${PREFIX}_metal_template.txt
```
```
module load R
R
require(dplyr)
require(data.table)

# Reformat to OR and P-value
i<-0
dfs = list()
file.names <- dir("/PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/METAL/", pattern=".*tbl$")
for(file in file.names) {
  i <- i+1
  data <- fread(file,header=T)
  beta <- data$Effect*(-1) # multiply by -1 if A1 isn't the minor allele 
  se <- data$StdErr
  P <- data$"P-value" %>% formatC(digits=4) %>% as.numeric()
  OR_all <- ifelse(is.na(P), "NA", paste(exp(beta) %>% round(digits=2), " (", exp(beta - 1.96*se) %>% round(digits=2), "-", exp(beta + 1.96*se) %>% round(digits=2), ")",sep=""))
  df <- data.frame(ORs = OR_all, Pval = P)
  prefix <- file %>% gsub(pattern="_metal1.tbl",replacement="")
  print(prefix)
  colnames(df) <- paste(prefix, colnames(df), sep = "_")
  df$SNP <- data$MarkerName
  dfs[[i]] <- df
}

# Merge all of the dataframes
library(tidyverse)
combined <- dfs %>% reduce(full_join, by = "SNP")
combined <- combined %>% relocate(SNP, .before=GS_ORs)
write.table(combined, file="regression_results.csv", quote=FALSE,row.names=F,sep=",")

# Get heterozygosity estimates here
het <- fread("NORMAL_metal1.tbl",header=T)
write.table(het, file="NORMAL_metal1.csv", quote=FALSE,row.names=F,sep=",")
q()
n

scp lakejs@biowulf.nih.gov://PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/METAL/regression_results.csv /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://PATH/TO/Julie/Julie_LRRK2_Condi/LRRK2_var_regression/METAL/NORMAL_metal1.csv /Users/lakejs/Desktop/
```
