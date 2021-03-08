## Testing the association of rs76904798 after removing LRRK2 coding variants

```
# Julie and Cornelis
# March 2021
```

### IPDGC
```
cd /data/LNG/Julie/Julie_LRRK2_Condi
mkdir no_coding_GWAS
cd no_coding_GWAS

# Recode the genotypes of interest as single allele dosage numbers 
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --snps 12:40629436, 12:40657700, 12:40671989, 12:40702911, 12:40707778, 12:40707861, 12:40713899, 12:40734202, 12:40740686 --recodeA --keep /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/NORMAL_covariates_GWAS.txt --out LRRK2_no_coding_GWAS
```

#### These are the variants we want to remove because they are the rare LRRK2 coding variants we analyzed
| MarkerName (hg19) | REF | ALT | RS-ID      | Region | Gene  | Protein Consequence |
|-------------------|-----|-----|------------|--------|-------|---------------------|
| 12:40629436       | T   | C   | rs33995463 | exonic | LRRK2 | L119P               |
| 12:40657700       | C   | G   | rs7308720  | exonic | LRRK2 | N551K               |
| 12:40671989       | A   | G   | rs10878307 | exonic | LRRK2 | I723V               |
| 12:40702911       | G   | A   | rs7133914  | exonic | LRRK2 | R1398H              |
| 12:40707778       | G   | A   | rs35507033 | exonic | LRRK2 | R1514Q              |
| 12:40707861       | C   | T   | rs33958906 | exonic | LRRK2 | P1542S              |
| 12:40713899       | T   | C   | rs35303786 | exonic | LRRK2 | M1646T              |
| 12:40734202       | G   | A   | rs34637584 | exonic | LRRK2 | G2019S              |
| 12:40740686       | A   | G   | rs33995883 | exonic | LRRK2 | N2081D              |

#### Create covariate files removing LRRK2 coding variants
```
module load R
R
require(dplyr)
require(data.table)
data <- read.table("LRRK2_no_coding_GWAS.raw",header=T)

# Remove all of the PD-linked rare LRRK2 coding variants
# Remove carriers of 551, 1398, 1646, G2019S, 2081
newdata1 <- subset(data, X12.40657700_G == 0 & X12.40702911_A == 0 & X12.40713899_C == 0 & X12.40734202_A == 0 & X12.40740686_G == 0)

# Remove all of the rare LRRK2 coding variants
coding_no_L119P <- subset(newdata1, X12.40671989_G == 0 & X12.40707778_A == 0 & X12.40707861_T == 0)

# Since L119P was not present in the MCGILL,MF,SPAIN3 and SPAIN4 cohorts due to low imputation quality, don't filter for this variant in those cohorts
# Add the dataset column to filter
cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T) %>% select(FID, DATASET)
`%notin%` <- Negate(`%in%`)
coding_no_L119P <- merge(coding_no_L119P, cov, by="FID")
yes_filter_L119P <- filter(coding_no_L119P, DATASET != "MCGILL" & DATASET != "MF" & DATASET != "SPAIN3" & DATASET != "SPAIN4" & DATASET != "GERMANY") %>% subset(X12.40629436_C == 0)
no_filter_L119P <- filter(coding_no_L119P, DATASET == "MCGILL" | DATASET == "MF" | DATASET == "SPAIN3" | DATASET == "SPAIN4" | DATASET == "GERMANY")
newdata2 <- rbind(yes_filter_L119P, no_filter_L119P)
newdata2$DATASET <- NULL

dim(newdata1) 
# 26579    15
dim(newdata2) 
# 20266    15

# Adding some additional sample info
cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T)
# Drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL
Mrg1 <- merge(newdata1,cov,by='FID')
Mrg2 <- merge(newdata2,cov,by='FID')
write.table(Mrg1, file="LRRK2_no_PDlinked_coding_with_COV.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Mrg2, file="LRRK2_no_coding_with_COV.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

### Make COV files for IPDGC
PCA_IPDGC() {
gwas_type=$1
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
mkdir $gwas_type
cd $gwas_type
cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do
    mkdir $line
    cd $line
    plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/${line}/${line} --keep /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/LRRK2_${gwas_type}_with_COV.txt --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
    plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
    plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
    plink --bfile prune --pca --out $line.LRRK2_condi_PCA_${gwas_type}
    # Send the .eigenvec files back to the working directory to combine into a new combined PC file
    scp $line.LRRK2_condi_PCA_${gwas_type}.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/${gwas_type}
    cd ..
    cat *${gwas_type}.eigenvec > ${gwas_type}_PCs.txt
done
}

PCA_IPDGC no_PDlinked_coding
PCA_IPDGC no_coding

### Create covariate files
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
R
require(data.table)
require(dplyr)
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS", pattern="LRRK2.*with_COV.txt") 
for(file in file.names) { 
  gwas_type <- file %>% gsub(pattern="LRRK2_", replacement="") %>% gsub(pattern="_with_COV.txt", replacement="")
  cov <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/", file,sep=""),header=T)
  PC <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/",gwas_type,"/",gwas_type,"_PCs.txt",sep=""),header=F)
  # We want to keep the phenotype information but get rid of the old PCs before merging
  cov2 <- cov[,c(1:21)]
  # Now add the new PCs
  Mrg <- merge(cov2,PC,by.x="FID",by.y="V1") 
  # Get rid of this column since it is a repeat of "FID"
  Mrg$V2 <- NULL 
  # Only keep the first 10 PCs
  Mrg2 <- Mrg[,c(1:31)]
  # Change the name of the first 10 PCs
  colnames(Mrg2)[22]  <- "PC1"
  colnames(Mrg2)[23]  <- "PC2"
  colnames(Mrg2)[24]  <- "PC3"
  colnames(Mrg2)[25]  <- "PC4"
  colnames(Mrg2)[26]  <- "PC5"
  colnames(Mrg2)[27]  <- "PC6"
  colnames(Mrg2)[28]  <- "PC7"
  colnames(Mrg2)[29]  <- "PC8"
  colnames(Mrg2)[30]  <- "PC9"
  colnames(Mrg2)[31]  <- "PC10"
  # The old sample selection files have an extra column, get rid of it 
  # Export the final coviariate files
  write.table(Mrg2, file=paste("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/",gwas_type,"/LRRK2_condi_covariates_",gwas_type,".txt",sep=""), quote=FALSE,row.names=F,sep="\t")
}
q()
n

subset_cohorts() {
gwas_type=$1
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/$gwas_type
cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt  | while read line
do
	# Pull out all of the lines that contain the cohort name
	# Combine these lines with the header
	grep -e FID -e $line LRRK2_condi_covariates_${gwas_type}.txt > LRRK2_condi_covariates_${gwas_type}.$line.txt
done
# Fix MF data...
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_${gwas_type}.MF.txt > temp
mv temp LRRK2_condi_covariates_${gwas_type}.MF.txt
# Fixed…
}
subset_cohorts no_PDlinked_coding
subset_cohorts no_coding
```

#### Run the GWAS
```
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 IPDGC_no_coding_GWAS.sh no_PDlinked_coding
sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 IPDGC_no_coding_GWAS.sh no_coding

# Reformat
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
reformat() {
gwas_type=$1
# Loop over to reformat the .hybrid files into .txt files
cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do
  cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/${gwas_type}/${gwas_type}_GWAS_CHR12/
  Rscript --vanilla /data/LNG/Julie/Julie_LRRK2_Condi/reformat_IPDGC.R ${gwas_type}_GWAS_CHR12.$line.PHENO_PLINK.glm.logistic.hybrid
done
}

reformat no_PDlinked_coding
reformat no_coding

# Organize the files 
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
organize() {
  gwas_type=$1
  cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/${gwas_type}/
  mkdir ${gwas_type}_COVARIATES
  mv *${gwas_type}.eigenvec ${gwas_type}_COVARIATES
  mv LRRK2_condi_covariates_${gwas_type}* ${gwas_type}_COVARIATES
  cd ${gwas_type}_GWAS_CHR12
  mkdir prep_files
  mv *.log prep_files
  mv *.hybrid prep_files
}

organize no_PDlinked_coding
organize no_coding
```

### UK Biobank

#### Create covariate files removing LRRK2 coding variants
```
# Determine LRRK2 carrier status
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
module load plink/2.3-alpha
plink2 --bgen /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb_imp_chr12_v3.bgen --snps rs33995463, rs7308720, rs10878307, rs7133914, rs35507033, rs33958906, rs35303786, rs34637584, rs33995883 --make-bed --sample /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb33601_imp_chr1_v3_s487395.sample --out LRRK2_no_coding_GWAS_UKB
module load plink
plink --bfile LRRK2_no_coding_GWAS_UKB --recodeA --out LRRK2_no_coding_GWAS_UKB2

cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
module load R
R
require(data.table)
require(dplyr)
PD_normal <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/UKB_PD_cases_control_over60.txt",header=T)
PD_normal <- PD_normal[,c(1:5)]
Proxy_normal <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/UKB_Proxy_cases_control_over60.txt",header=T)
Proxy_normal <- Proxy_normal[,c(1:5)]
LRRK2_status <- fread("LRRK2_no_coding_GWAS_UKB2.raw",header=T)
LRRK2_status$IID <- NULL
LRRK2_status$PAT <- NULL
LRRK2_status$MAT <- NULL
LRRK2_status$SEX <- NULL
LRRK2_status$PHENOTYPE <- NULL

# Add LRRK2 status to the final PD and Proxy datasets
PD_normal_LRRK2 <- merge(PD_normal,LRRK2_status,by.x='FID',by.y='FID')
Proxy_normal_LRRK2 <- merge(Proxy_normal,LRRK2_status,by.x='FID',by.y='FID')

PD_no_PDlinked_coding <- subset(PD_normal_LRRK2, rs7308720_G == 0 & rs7133914_A == 0 & rs35303786_C == 0 & rs34637584_A == 0 & rs33995883_G == 0)
PD_no_coding <- subset(PD_no_PDlinked_coding, rs33995463_C == 0 & rs10878307_G == 0 & rs35507033_A == 0 & rs33958906_T == 0)
Proxy_no_PDlinked_coding <- subset(Proxy_normal_LRRK2, rs7308720_G == 0 & rs7133914_A == 0 & rs35303786_C == 0 & rs34637584_A == 0 & rs33995883_G == 0)
Proxy_no_coding <- subset(Proxy_no_PDlinked_coding, rs33995463_C == 0 & rs10878307_G == 0 & rs35507033_A == 0 & rs33958906_T == 0)

write.table(PD_no_PDlinked_coding, file="UKB_PD_cases_control_over60_no_PDlinked_coding.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_no_coding, file="UKB_PD_cases_control_over60_no_coding.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Proxy_no_PDlinked_coding, file="UKB_Proxy_cases_control_over60_no_PDlinked_coding.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Proxy_no_coding, file="UKB_Proxy_cases_control_over60_no_coding.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

# Calculate PCs
module load flashpca
module load plink

cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS
mkdir UKB_no_coding_GWAS
cd UKB_no_coding_GWAS

sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/flashpca_UKB_no_coding.sh no_PDlinked_coding
sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/flashpca_UKB_no_coding.sh no_coding

# Merge files in R
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS
module load R
R
require(data.table)
require(dplyr)
# Load the full covariates file
cov <- fread("/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt",header=T)
# Remove the PC columns from cov so that we can merge with the new PCs
cov2 <- cov %>% select(1:8)
# Pull the subset phenotype files from your working directory
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS", pattern="^UKB_P")
for(file in file.names) {  
  pc <- fread(paste("pcs_",file,sep=""),header=T)
  pc$IID <- NULL
  pheno <- fread(paste("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/",file,sep=""),header=T)
  # This info will become redundant when merging, keep (PD/CONTROL/PROXY) STATUS and LRRK2 carrier status
  pheno$IID <- pheno$SEX <- pheno$AGE <- NULL
  Mrg = merge(cov2,pheno,by='FID')
  Mrg2 = merge(Mrg,pc,by='FID')
  write.table(Mrg2, file=paste("COV_",file,sep=""), quote=FALSE,row.names=F,sep="\t")
}
q()
n

# Replace PD/CONTROL/PROXY STATUS with numbers
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS
ls COV_UKB_P* > UKB_COV_files.txt

cat UKB_COV_files.txt | while read line
do 
  sed -i 's/PD/2/g' $line
  sed -i 's/PROXY/2/g' $line
  sed -i 's/CONTROL/1/g' $line
done

# Organize and remove some files
mkdir flashpca_files
mv eigenvalues_UKB_P* flashpca_files
mv eigenvectors_UKB_P* flashpca_files
mv pcs_UKB_P* flashpca_files
mv pve_UKB_P* flashpca_files
rm *bed
rm *bim
rm *fam
rm *log
rm pruned_data*
```

#### Perform GWAS
```
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS
module load plink/2.0-dev-20191128
sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 UKB_GWAS_no_coding.sh no_PDlinked_coding
sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 UKB_GWAS_no_coding.sh no_coding

# Reformat plink2 GWAS output
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS
ls COV_UKB_P*.hybrid > GWAS_files.txt
module load python/3.6
cat GWAS_files.txt | while read line
do 
  # Filter the results by A1_FREQ (minor allele frequency) >=0.0001 --> to input.txt \
  awk '{ if($13 >= 0.0001) { print }}' $line > input.txt
  # Reformat the plink2 results for meta-analysis using python
  if [[ $line == "COV_UKB_PD"* ]]; then
    python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
    --outfile toMeta.${line%%.*}.txt --B-or-C B
  elif [[ $line == "COV_UKB_Proxy"* ]]; then
    python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
    --outfile toProxy.${line%%.*}.txt --B-or-C B; fi
done

###Adjust Proxy cases to the same scale as PD cases
# Make the toProxy files into .csv files
module load R
R
require(data.table)
require(dplyr)
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS", pattern="toProxy") 
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
    python /data/CARD/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
    --infile $line --beta-proxy beta --se-proxy se --p-proxy P --outfile $newname
done

# Convert the "normalized" proxy .csv files back to .txt files
R
require(data.table)
require(dplyr)
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS", pattern="toMeta.COV_UKB_Proxy.*csv") 
for(file in file.names) {  
  data <- fread(file,header=T)
  newname <- file %>% gsub(pattern=".csv", replacement=".txt")
  write.table(data, file=newname,quote=F,row.names=F,sep="\t")
}
q()
n

# Reformat and move to the same directory as the IPDGC CHR12 files for meta-analysis
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/UKB_no_coding_GWAS

ls toMeta*.txt > meta_files.txt
cat meta_files.txt | while read line
do
    # Get rid of the :N:N off of the markername so that it matches IPDGC
    sed -i -r 's/:[ACGT]+//g' $line
    pattern1="*_over60_"
    pattern2=".txt"
    replacement=""
    newname1="${line/$pattern1/$replacement}"
    newname2="${newname1/$pattern2/$replacement}"
    if [[ $line == "toMeta.COV_UKB_PD"* ]]; then
        cut -f1-7 $line | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/${newname2}/${newname2}_GWAS_CHR12/${newname2}_GWAS_CHR12.UKBPD.txt
    elif [[ $line == "toMeta.COV_UKB_Proxy"* ]]; then
        cut -f1,2,3,7,15,16,17 $line | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/${newname2}/${newname2}_GWAS_CHR12/${newname2}_GWAS_CHR12.UKBproxy.txt; fi
done
```
### Forest plot of rs76904798 association after removing PD-linked and all rare LRRK2 coding variants 
```
cd /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS

head -1 /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.DUTCH.txt > header.txt
grep "12:40614434" /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL*txt > temp
cat header.txt temp > NORMAL_12:40614434.txt
sed -e 's/.*GWAS_CHR12.//g' NORMAL_12:40614434.txt |  sed -e 's/.txt:12:40614434//g' > NORMAL_12:40614434v2.txt

grep "12:40614434" /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/no_coding/no_coding_GWAS_CHR12/no_coding_GWAS_CHR12.*.txt > temp
cat header.txt temp > no_coding_12:40614434.txt
sed -e 's/.*GWAS_CHR12.//g' no_coding_12:40614434.txt | sed -e 's/.txt:12:40614434//g' > no_coding_12:40614434v2.txt

grep "12:40614434" /data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/no_PDlinked_coding/no_PDlinked_coding_GWAS_CHR12/no_PDlinked_coding_GWAS_CHR12.*.txt > temp
cat header.txt temp > no_PDlinked_coding_12:40614434.txt
sed -e 's/.*GWAS_CHR12.//g' no_PDlinked_coding_12:40614434.txt |  sed -e 's/.txt:12:40614434//g' > no_PDlinked_coding_12:40614434v2.txt

R
require(dplyr)
library(metafor)
data_normal <- read.table("NORMAL_12:40614434v2.txt", header = T)
data_no_coding <- read.table("no_coding_12:40614434v2.txt", header = T)
data_no_PDlinked <- read.table("no_PDlinked_coding_12:40614434v2.txt", header = T)

labs_normal <- gsub(".*\\.","", data_normal$ID)
labs_normal <- gsub("NEUROX_DBGAP", "NEUROX", labs_normal)
labs_normal <- gsub("UKBPD", "UKBIO_case", labs_normal)
labs_normal <- gsub("UKBproxy", "UKBIO_proxy", labs_normal)
yi_normal   <- data_normal$beta
sei_normal  <- data_normal$LOG.OR._SE
resFe_normal  <- rma(yi=yi_normal, sei=sei_normal, method="FE")
resRe_normal  <- rma(yi=yi_normal, sei=sei_normal)

yi_no_coding <- data_no_coding$beta
sei_no_coding <- data_no_coding$LOG.OR._SE
resFe_no_coding <- rma(yi=yi_no_coding, sei=sei_no_coding, method="FE")
resRe_no_coding<- rma(yi=yi_no_coding, sei=sei_no_coding)

yi_no_PDlinked <- data_no_PDlinked$beta
sei_no_PDlinked <- data_no_PDlinked$LOG.OR._SE
resFe_no_PDlinked <- rma(yi=yi_no_PDlinked, sei=sei_no_PDlinked, method="FE")
resRe_no_PDlinked  <- rma(yi=yi_no_PDlinked, sei=sei_no_PDlinked)

pdf(file = "12:40614434_combined.pdf", width = 8, height = 7)
Pvalue_normal <- formatC(resFe_normal$pval, digits=4)
Pvalue_no_coding <- formatC(resFe_no_coding$pval, digits=4)
Pvalue_no_PDlinked <- formatC(resFe_no_PDlinked$pval, digits=4)

#Make it so that all datasets are included even if NA for the variant
options(na.action = "na.pass")

par(mfrow=c(1,3))

par(mar=c(5,4,1,1))
forest(resFe_normal, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""),
       slab=labs_normal, mlab="Fixed Effects", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1, "Unconditioned", cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_normal,sep=""), cex=1.2, font=2)

par(mar=c(5,0,1,1))
forest(resFe_no_PDlinked, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=NA, mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1, expression(bold(Delta)~bolditalic("LRRK2 ")~bold("PD-linked")), cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_no_PDlinked,sep=""), cex=1.2, font=2)
#adding this for the title
text(0, 18, "rs76904798", cex=1.5, font=2)

par(mar=c(5,0,1,2))
forest(resFe_no_coding, annotate=TRUE, xlim=c(-2.25,3.25), width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=NA, mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1,expression(bold(Delta)~bolditalic("LRRK2 ")~bold("Rare")), cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_no_coding,sep=""), cex=1.2, font=2)
dev.off()
q()
n

scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/no_coding_GWAS/12:40614434_combined.pdf /Users/lakejs/Desktop/
```

Main conclusion => rs76904798 association signal is independent of LRRK2 coding variants
