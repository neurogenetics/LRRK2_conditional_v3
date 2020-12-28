# LRRK2 conditional 2020

`LNG ❤️ Open Science 😍`

 - **Project:** LRRK2 conditional GWAS
 - **Author(s):** Cornelis Blauwendraat, Julie Lake, Hampton Leonard (LNG) Nicole Bryant-Garner (Duke Uni)
 - **Date Last Updated:** November 2020
    - **Update Description:** Edits to README

---
### Quick Description: 
LRRK2 is an important gene for PD and both rare highly damaging missense variants (including G2019S) and common non-coding variants (rs76904798) have been associated with PD risk. Rare highly damaging missense variants often results in a very high increase of risk for PD >10 OR, while the common non-coding variant have the usual moderate risk increase of ~1.15. The goal here is to identify whether other more common coding variants also result in an higher risk of PD or whether they are just associated with PD as result of complex linkage disequilibrium (LD) structures.


### Motivation/Goals:
1) Understanding the underlying data and creating an overview of the data...
2) Perform quick and dirty GWAS excluding risk and G2019S variant...
3) Perform GWAS excluding risk and G2019S variant on a cohort level basis and meta-analyze

### Link to Manuscript:
TBD

## Structure of README:
### [1. Understanding the underlying data and creating an overview of the data](#1-Understanding-the-underlying-data-and-creating-an-overview-of-the-data)
This section goes through:
- Intro to the IPDGC data
- Checking the imputation quality of the IPDGC data
- Assessing frequency of LRRK2 G2019S and rs76904798 in the data
- Overview of the full data and selection of data

### [2. Create covariate files for each IPDGC cohort](#2-Create-covariate-files-for-each-IPDGC-cohort)
This section goes through:
- Creating new PC's for each cohort
- Creating covariate files

### [3. Perform cohort-level GWAS on IPDGC data excluding rs76904798 and G2019S](#3-Perform-cohort-level-GWAS-on-IPDGC-data-excluding-rs76904798-and-G2019S)
This section goes through: 
- Extracting all LRRK2 coding variants
- Performing GWAS of chromosome 12 for each cohort
- Prepping before meta-analysis

### [4. Adding in UKBiobank](#4-Adding-in-UKBiobank)
 This section goes through: 
- Subsetting the UK Biobank data
- Making covariate files
- Performing GWAS on CHR12
- Reformatting data for meta-analysis

### [5. Making forest plots for LRRK2 coding variants](#5-Making-forest-plots-for-LRRK2-coding-variants)
This section goes through: 
- Pulling the amino acid changes for each variant
- Using metafor to make forest plots for LRRK2 coding variants

### [6. Check LD co-inheritance of LRRK2 coding variants](#6-Check-LD-co-inheritance-of-LRRK2-coding-variants)
This section goes through:
- Making frequency files for IPDGC and UKB data
- Checking co-inheritance of G2019S, rs76904798 and N2081D with all other coding variants
- Determining a frequency cutoff for LRRK2 coding variants

### [7. Make final tables and figures](#7-Make-final-tables-and-figures)
This section goes through:
- Preparing tables for manuscript
- Preparing figures for manuscript

---
## 1. Understanding the underlying data and creating an overview of the data
This section goes through:
- Intro to the data
- Assessing frequency of LRRK2 G2019S and rs76904798 in the data
- Checking the imputation quality of the data
- Overview of the full data and selection of which data to continue with

### 1.1 - Intro to the IPDGC data

Path to working folder: /data/LNG/Julie/Julie_LRRK2_Condi

Path to IPDGC genetics data: /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins.bim/bed/fam

This data is filtered for a lot of things...check https://github.com/neurogenetics/GWAS-pipeline for more details plus variants are filtered for a very conservative R2 > 0.8 and data is filtered for relatedness in the full dataset for pihat < 0.125

```
# Variants of interest:
# LRRK2 G2019S => hg19 12:40734202:G:A
# rs76904798 => hg19 12:40614434:C:T

cd /data/LNG/Julie/Julie_LRRK2_Condi 
module load plink

# Simple test
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202,12:40614434 --assoc --out test
# Among remaining phenotypes, 21478 are cases and 24388 are controls.
 CHR           SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  12   12:40614434   40614434    T   0.1564   0.1406    C        45.29    1.702e-11        1.133 
  12   12:40734202   40734202    A 0.007297 0.0005648    G        203.6    3.346e-46        13.01 
# Confirming associations, so thats good :)

rm test* 
```

### 1.2 - Checking the imputation quality of the IPDGC data

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

# Make cohort loop over file
echo NEUROX_DBGAP > cohort_file.txt
echo MCGILL >> cohort_file.txt
echo VANCE >> cohort_file.txt
echo NIA >> cohort_file.txt
echo GERMANY >> cohort_file.txt
echo PPMI >> cohort_file.txt
echo SPAIN3 >> cohort_file.txt
echo HBS >> cohort_file.txt
echo SHULMAN >> cohort_file.txt
echo MF >> cohort_file.txt
echo DUTCH >> cohort_file.txt
echo PDBP >> cohort_file.txt
echo SPAIN4 >> cohort_file.txt

# Pull the imputation quality information for each cohort for variants of interest
cat cohort_file.txt | while read line
do 
	zless /data/LNG/CORNELIS_TEMP/PD_AAO/${line}/chr12.info.gz | grep 12:40614434 | cut -f7,8 >> temp1.txt
	zless /data/LNG/CORNELIS_TEMP/PD_AAO/${line}/chr12.info.gz | grep 12:40734202 | cut -f7,8 >> temp2.txt
	zless /data/LNG/CORNELIS_TEMP/PD_AAO/${line}/chr12.info.gz | grep 12:123326598 | cut -f7,8 >> temp3.txt
	zless /data/LNG/CORNELIS_TEMP/PD_AAO/${line}/chr12.info.gz | grep 12:40740686 | cut -f7,8 >> temp4.txt
done

module load R
R
require(data.table)
require(dplyr)

data1 <- fread("temp1.txt",header=F)
data2 <- fread("temp2.txt",header=F)
data3 <- fread("temp3.txt",header=F)
data4 <- fread("temp4.txt",header=F)

colnames(data1) <- c("12:40614434","How")
colnames(data2) <- c("12:40734202","How")
colnames(data3) <- c("12:123326598","How")
colnames(data4) <- c("12:40740686","How")

cohorts <- fread("cohort_file.txt",header=F)
colnames(cohorts) <- c("Dataset")

data <- cbind(cohorts,data1,data2,data3,data4)
write.table(data,file="Imputation_quality.txt",quote=FALSE,row.names=F,sep="\t")
q()
n

rm temp*.txt

# Copy file home
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/Imputation_quality.txt /Users/lakejs/Desktop/

# Summary
# rs76904798 => hg19 12:40614434:C:T -> very well imputed, present in almost all data
# LRRK2 G2019S => hg19 12:40734202:G:A -> also not bad...
# rs10847864 => 12:123326598:G:T (positive control) -> some bad cohorts 
# LRRK2 N2081D => 12:40740686:A:G -> very well imputed, present in almost all data
```

#### IPDGC imputation quality table:

| Dataset      | 12:40614434 | How       | 12:40734202 | How       | 12:123326598 | How       | 12:40740686 | How       |
|--------------|-------------|-----------|-------------|-----------|--------------|-----------|-------------|-----------|
| NEUROX_DBGAP | 0.99998     | Genotyped | 0.98795     | Genotyped | 0.99841      | Genotyped | 0.9977      | Genotyped |
| MCGILL       | 1           | Genotyped | 0.91193     | Imputed   | 0.99912      | Genotyped | 0.94829     | Imputed   |
| VANCE        | 1           | Genotyped | 0.99384     | Imputed   | 0.90858      | Imputed   | 0.9995      | Imputed   |
| NIA          | 0.99775     | Imputed   | 0.93746     | Imputed   | 0.87344      | Imputed   | 0.9903      | Imputed   |
| GERMANY      | 0.9961      | Imputed   | 0.85943     | Imputed   | 0.87789      | Imputed   | 0.99355     | Imputed   |
| PPMI         | 0.99999     | Genotyped | 0.97352     | Genotyped | 0.99941      | Genotyped | 0.99474     | Genotyped |
| SPAIN3       | 0.98578     | Imputed   | 0.99945     | Genotyped | **0.76704**      | Imputed   | 0.99982     | Genotyped |
| HBS          | 0.99999     | Genotyped | 0.99982     | Genotyped | 0.99935      | Genotyped | 0.99999     | Genotyped |
| SHULMAN      | 0.99557     | Genotyped | 0.99997     | Genotyped | 0.84249      | Imputed   | 0.99989     | Genotyped |
| MF           | 0.99684     | Imputed   | 0.9629      | Imputed   | **0.77759**      | Imputed   | 0.98894     | Imputed   |
| DUTCH        | 0.99903     | Imputed   | 0.87355     | Imputed   | 0.90051      | Imputed   | 0.99208     | Imputed   |
| PDBP         | 0.99999     | Genotyped | 0.97688     | Genotyped | 0.99906      | Genotyped | 0.99999     | Genotyped |
| SPAIN4       | 0.98652     | Imputed   | 0.99914     | Genotyped | **0.76579**      | Imputed   | 0.99936     | Genotyped |


#### Add the positive control rs10847864 back to the hardcalls

```
# See that in MF, SPAIN3 and SPAIN4, rs10847864 has an imputation quality below the hardcall threshold of 0.8 
# Manually add this variant back to the hardcalls

# Make plink binary files for the variant in each of MF, SPAIN3, SPAIN4 datasets
cd /data/LNG/Julie/Julie_LRRK2_Condi
module load plink

# Pull the variant from the cohorts where it's missing due to low imputation quality
for cohort in {"MF","SPAIN3","SPAIN4"};
do
	plink --vcf /data/LNG/CORNELIS_TEMP/PD_AAO/${cohort}/chr12.dose.vcf.gz --make-bed --out s1 --double-id
	plink --bfile s1 --snps 12:123326598 --make-bed --out ${cohort}_rs10847864_only
done

# Merge these new files with the HARDCALLS 
# Note this will take a while
sinteractive --mem=240g --cpus-per-task=20
module load plink
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --bmerge MF_rs10847864_only --out hardcalls1
plink --bfile hardcalls1 --bmerge SPAIN3_rs10847864_only --out hardcalls2
plink --bfile hardcalls2 --bmerge SPAIN4_rs10847864_only --out HARDCALLS_with_rs10847864

rm hardcalls1*
rm hardcalls2*
mkdir rs10847864_prep
mv SPAIN3* rs10847864_prep
mv SPAIN4* rs10847864_prep
mv MF* rs10847864_prep

# Use HARDCALLS_with_rs10847864.bed/bim/fam for future analysis
```

### 1.3 - Assessing frequency of LRRK2 G2019S, rs76904798 and N2081D in the data
```
# Check allelic distribution
cd /data/LNG/Julie/Julie_LRRK2_Condi
plink --bfile HARDCALLS_with_rs10847864 --snps 12:40734202,12:40614434,12:40740686 --model --out allelic_dist

R
require(dplyr)
require(data.table)
data <- fread("allelic_dist.model",header=T)
data2 <- subset(data, TEST=="GENO")
write.table(data2, file="allelic_dist.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

cat allelic_dist.txt
CHR	SNP	A1	A2	TEST	AFF	UNAFF	CHISQ	DF	P
12	12:40614434	T	C	GENO	558/5602/15318	486/5885/18017	46.02	2	1.014e-10
12	12:40734202	A	G	GENO	1/236/16072	0/20/17685	NA	NA	NA
12	12:40740686	G	A	GENO	11/992/20475	13/896/23479	25.83	2	2.46e-06
```

### 1.4 - Overview of the full data and selection of data

```
# Recode the genotypes of interest as single allele dosage numbers 

# Use for the no rs76904798 + no G2019S dataset
plink --bfile HARDCALLS_with_rs10847864 --snps 12:40734202,12:40614434 --recodeA --out LRRK2_condi_variant_selection

# Use for the no N2081D + no G2019S dataset
# Call this the "special" conditional GWAS: see if rs76904798 signal remains after removal of N2081D
plink --bfile HARDCALLS_with_rs10847864 --snps 12:40734202,12:40740686 --recodeA --out LRRK2_condi_special_variant_selection
```

#### Subset the data to include only homozygous reference carriers of both variants

```
module load R
R
require(dplyr)
require(data.table)
data <- read.table("LRRK2_condi_variant_selection.raw",header=T)
data2 <- read.table("LRRK2_condi_special_variant_selection.raw",header=T)

# X12.40614434_T is rs76904798 and X12.40734202_A is G2019S
newdata <- subset(data, X12.40614434_T == 0 & X12.40734202_A == 0) 
dim(newdata) 
# 24532     8

# X12:40740686_G is N2081D and X12.40734202_A is G2019S
newdata2 <- subset(data2, X12.40740686_G == 0 & X12.40734202_A == 0)
dim(newdata2) 
# 32227     8

# Add some additional sample info
cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T)

# Drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL

Mrg = merge(newdata,cov,by='FID')
dim(Mrg) 
# 24532    44
Mrg2 = merge(newdata2,cov,by='FID')
dim(Mrg2) 
# 32227    44

# Datasets with good data (the ones with homo-ref carriers):
Mrg$DATASET %>% unique()
[1] "NEUROX_DBGAP" "MCGILL"       "VANCE"        "NIA"          "GERMANY"     
[6] "PPMI"         "SPAIN3"       "HBS"          "SHULMAN"      "MF"          
[11] "DUTCH"        "PDBP"         "SPAIN4"      

Mrg2$DATASET %>% unique()
[1] "NEUROX_DBGAP" "MCGILL"       "VANCE"        "NIA"          "GERMANY"     
[6] "PPMI"         "SPAIN3"       "HBS"          "SHULMAN"      "MF"          
[11] "DUTCH"        "PDBP"         "SPAIN4"    

# Create a file for selecting individuals who are homo-ref carriers
write.table(Mrg, file="LRRK2_condi_sample_selection.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Mrg2, file="LRRK2_condi_special_sample_selection.txt", quote=FALSE,row.names=F,sep="\t")

# Display the case-control distribution for each dataset
Mrg_grouped <- Mrg %>% group_by(DATASET) %>% summarise(Case = sum(PHENOTYPE == 2), Control = sum(PHENOTYPE == 1), TOTAL = n()) %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "SUM"))
write.table(Mrg_grouped, file="LRRK2_condi_sample_selection_grouped.txt", quote=FALSE,row.names=F,sep="\t")

Mrg2_grouped <- Mrg2 %>% group_by(DATASET) %>% summarise(Case = sum(PHENOTYPE == 2), Control = sum(PHENOTYPE == 1), TOTAL = n()) %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "SUM"))
write.table(Mrg2_grouped, file="LRRK2_condi_special_sample_selection_grouped.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

# Copy files
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_sample_selection_grouped.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_special_sample_selection_grouped.txt /Users/lakejs/Desktop/
```

<table>
<tr><th>Data for conditional GWAS (no rs76904798 + no G2019S)</th><th>Data for special conditional GWAS (no N2081D + no G2019S)</th></tr>
<tr><td>

| DATASET      | Case  | Control | TOTAL |
|--------------|-------|---------|-------|
| DUTCH        | 539   | 1476    | 2015  |
| GERMANY      | 505   | 723     | 1228  |
| HBS          | 365   | 333     | 698   |
| MCGILL       | 410   | 669     | 1079  |
| MF           | 537   | 565     | 1102  |
| NEUROX_DBGAP | 3801  | 4248    | 8049  |
| NIA          | 579   | 2233    | 2812  |
| PDBP         | 384   | 204     | 588   |
| PPMI         | 252   | 134     | 386   |
| SHULMAN      | 553   | 152     | 705   |
| SPAIN3       | 1459  | 989     | 2448  |
| SPAIN4       | 1634  | 1147    | 2781  |
| VANCE        | 423   | 218     | 641   |
| SUM          | 11441 | 13091   | **24532** |

</td><td>

| DATASET      | Case  | Control | TOTAL |
|--------------|-------|---------|-------|
| DUTCH        | 719   | 1905    | 2624  |
| GERMANY      | 668   | 902     | 1570  |
| HBS          | 489   | 447     | 936   |
| MCGILL       | 554   | 890     | 1444  |
| MF           | 713   | 739     | 1452  |
| NEUROX_DBGAP | 5036  | 5500    | 10536 |
| NIA          | 785   | 2898    | 3683  |
| PDBP         | 486   | 270     | 756   |
| PPMI         | 339   | 162     | 501   |
| SHULMAN      | 719   | 187     | 906   |
| SPAIN3       | 1969  | 1288    | 3257  |
| SPAIN4       | 2219  | 1490    | 3709  |
| VANCE        | 569   | 284     | 853   |
| SUM          | 15265 | 16962   | **32227** |

</td></tr></table>

## 2. Create covariate files for each IPDGC cohort

This section goes through:
- Creating new PC's for each cohort
- Creating covariate files
	- Total cohorts to include is 13 
	- For each cohort, make covariate files for: 
		Normal GWAS: PD vs control 
		Conditional GWAS: PD vs control (no rs76904798 + no G2019S)
		Special conditional GWAS: PD vs control (no N2081D + no G2019S)
	- Total covariate files needed is 13 x 3 = 39

Also note that most data was already preprocessed according to this https://github.com/neurogenetics/GWAS-pipeline 
and has been used in previous GWAS such as https://pubmed.ncbi.nlm.nih.gov/31701892/ and https://pubmed.ncbi.nlm.nih.gov/30957308/

### 2.1 Create new PC's for each cohort

First loop over all cleaned unimputed data to create new PC's

#### PC files without G2019S and rs76904798

```
### Loop for making PCA
cd /data/LNG/Julie/Julie_LRRK2_Condi
cat cohort_file.txt | while read line
do 
	# First copy over the pre_impute_vcf_files to new directories
	mkdir $line
	cp /data/LNG/CORNELIS_TEMP/PD_AAO/pre_impute_vcf_files/$line/{$line.bed,$line.bim,$line.fam} /data/LNG/Julie/Julie_LRRK2_Condi/$line
	
	cd $line
	plink --bfile $line --keep /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_sample_selection.txt --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_CONDI
	
	# Send the .eigenvec files back to the working directory to combine into a new combined PC file	
	scp $line.LRRK2_condi_PCA_CONDI.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	
	cd ..
done
```

#### PC files without G2019S and N2081D

```
### Loop for making PCA
cat cohort_file.txt | while read line
do 
	cd $line
	plink --bfile $line --keep /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_special_sample_selection.txt --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_SPECIAL
	
	# Send the .eigenvec files back to the working directory to combine into a new combined PC file
	scp $line.LRRK2_condi_PCA_SPECIAL.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	cd ..
done
```

#### PC files with G2019S, N2081D and rs76904798 <= meaning normal files...

```
## Loop for making PCA
# All samples!
cat cohort_file.txt  | while read line
do 
	cd $line
	# The only difference here is that we don't use --keep to filter out people with the variants
	plink --bfile $line --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_NORMAL
	
	# Send the .eigenvec files back to the working directory to combine into a new combined PC file
	scp $line.LRRK2_condi_PCA_NORMAL.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	cd ..
done
```

### 2.2 Create covariate files

#### Merge new PC's in R with other phenotype data

```
# Combine the eigenvec files for each GWAS
cat *NORMAL.eigenvec > NORMAL_PCs.txt
cat *SPECIAL.eigenvec > SPECIAL_PCs.txt
cat *CONDI.eigenvec > CONDI_PCs.txt

module load R
R
NORMAL_cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T)
SPECIAL_cov <- read.table("LRRK2_condi_special_sample_selection.txt",header=T)
CONDI_cov <- read.table("LRRK2_condi_sample_selection.txt",header=T)

NORMAL_PC <- read.table("NORMAL_PCs.txt",header=F)
SPECIAL_PC <- read.table("SPECIAL_PCs.txt",header=F)
CONDI_PC <- read.table("CONDI_PCs.txt",header=F)

# We want to keep the phenotype information but get rid of the old PCs before merging
NORMAL_cov2 <- NORMAL_cov[,c(1:10)]
SPECIAL_cov2 <- SPECIAL_cov[,c(1:14)]
CONDI_cov2 <- CONDI_cov[,c(1:14)]

# Now add the new PCs
NORMAL_Mrg <- merge(NORMAL_cov2,NORMAL_PC,by.x="FID",by.y="V1") 
SPECIAL_Mrg <- merge(SPECIAL_cov2,SPECIAL_PC,by.x="FID",by.y="V1") 
CONDI_Mrg <- merge(CONDI_cov2,CONDI_PC,by.x="FID",by.y="V1") 

# Get rid of this column since it is a repeat of "FID"
NORMAL_Mrg$V2 <- NULL 
SPECIAL_Mrg$V2 <- NULL 
CONDI_Mrg$V2 <- NULL 

# Only keep the first 10 PCs
NORMAL_Mrg2 <- NORMAL_Mrg[,c(1:20)]
SPECIAL_Mrg2 <- SPECIAL_Mrg[,c(1:24)]
CONDI_Mrg2 <- CONDI_Mrg[,c(1:24)]

# Change the name of the first 10 PCs

# Do this for the normal files
colnames(NORMAL_Mrg2)[11]  <- "PC1"
colnames(NORMAL_Mrg2)[12]  <- "PC2"
colnames(NORMAL_Mrg2)[13]  <- "PC3"
colnames(NORMAL_Mrg2)[14]  <- "PC4"
colnames(NORMAL_Mrg2)[15]  <- "PC5"
colnames(NORMAL_Mrg2)[16]  <- "PC6"
colnames(NORMAL_Mrg2)[17]  <- "PC7"
colnames(NORMAL_Mrg2)[18]  <- "PC8"
colnames(NORMAL_Mrg2)[19]  <- "PC9"
colnames(NORMAL_Mrg2)[20]  <- "PC10"

# Do this for the special conditional files
colnames(SPECIAL_Mrg2)[15]  <- "PC1"
colnames(SPECIAL_Mrg2)[16]  <- "PC2"
colnames(SPECIAL_Mrg2)[17]  <- "PC3"
colnames(SPECIAL_Mrg2)[18]  <- "PC4"
colnames(SPECIAL_Mrg2)[19]  <- "PC5"
colnames(SPECIAL_Mrg2)[20]  <- "PC6"
colnames(SPECIAL_Mrg2)[21]  <- "PC7"
colnames(SPECIAL_Mrg2)[22]  <- "PC8"
colnames(SPECIAL_Mrg2)[23]  <- "PC9"
colnames(SPECIAL_Mrg2)[24]  <- "PC10"

# Do this for the conditional files
colnames(CONDI_Mrg2)[15]  <- "PC1"
colnames(CONDI_Mrg2)[16]  <- "PC2"
colnames(CONDI_Mrg2)[17]  <- "PC3"
colnames(CONDI_Mrg2)[18]  <- "PC4"
colnames(CONDI_Mrg2)[19]  <- "PC5"
colnames(CONDI_Mrg2)[20]  <- "PC6"
colnames(CONDI_Mrg2)[21]  <- "PC7"
colnames(CONDI_Mrg2)[22]  <- "PC8"
colnames(CONDI_Mrg2)[23]  <- "PC9"
colnames(CONDI_Mrg2)[24]  <- "PC10"

# Export the final coviariate files
write.table(NORMAL_Mrg2, file="LRRK2_condi_covariates_NORMAL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(SPECIAL_Mrg2, file="LRRK2_condi_covariates_SPECIAL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(CONDI_Mrg2, file="LRRK2_condi_covariates_CONDI.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Subset each cohort due to potentially overlapping sample names

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

cat cohort_file.txt  | while read line
do
	# Pull out all of the lines that contain the cohort name
	# Combine these lines with the header
	grep -e FID -e $line LRRK2_condi_covariates_NORMAL.txt > LRRK2_condi_covariates_NORMAL.$line.txt
	grep -e FID -e $line LRRK2_condi_covariates_SPECIAL.txt > LRRK2_condi_covariates_SPECIAL.$line.txt
	grep -e FID -e $line LRRK2_condi_covariates_CONDI.txt > LRRK2_condi_covariates_CONDI.$line.txt
done


# Fix MF data...
# The MF covariates file also pulled out other cohorts since "MF" is contained within other cohorts' IDs (e.g. HBS_PD_INVDG562MF3)

# Want to return all non matching lines -- lines without HBS, PDBP, SPAIN4
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_NORMAL.MF.txt > temp
mv temp LRRK2_condi_covariates_NORMAL.MF.txt
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_SPECIAL.MF.txt > temp
mv temp LRRK2_condi_covariates_SPECIAL.MF.txt
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_CONDI.MF.txt > temp
mv temp LRRK2_condi_covariates_CONDI.MF.txt

# Fixed…
```

#### Reorganize the files

```
# Move the NORMAL covariate files into a new directory
mkdir NORMAL_COVARIATES
mv *NORMAL.eigenvec NORMAL_COVARIATES
mv LRRK2_condi_covariates_NORMAL* NORMAL_COVARIATES

# Move the SPECIAL covariate files into a new directory
mkdir SPECIAL_COVARIATES
mv *SPECIAL.eigenvec SPECIAL_COVARIATES
mv LRRK2_condi_covariates_SPECIAL* SPECIAL_COVARIATES

# Move the conditional covariate files into a new directory
mkdir CONDI_COVARIATES
mv *CONDI.eigenvec CONDI_COVARIATES
mv LRRK2_condi_covariates_CONDI* CONDI_COVARIATES
```

## 3. Perform cohort-level GWAS on IPDGC data excluding rs76904798 and G2019S

This section goes through: 
- Extracting all LRRK2 coding variants
- Performing GWAS of chromosome 12 for each cohort
- Prepping before meta-analysis

### 3.1 Extract all LRRK2 coding variants

```
cd /data/LNG/Julie/Julie_LRRK2_Condi
mkdir HRC_LRRK2
cd HRC_LRRK2

head -1 /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > header.txt

# Pull out the variants from CHR12, then filter by LRRK2, exonic, nonsynonymous --> LRRK2 coding variants
awk '$1 == "12"' /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt | grep LRRK2 | grep exonic | sed 's/ /_/g' | grep nonsynonymous > temp
cat header.txt temp > LRRK2_HRC_coding.txt

# Manually add in rs76904798 (non-coding)
grep rs76904798 /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > risk_variant.txt
cat LRRK2_HRC_coding.txt risk_variant.txt > LRRK2_HRC_coding_V2.txt

# Manually add in non-LRRK2 CHR12 GWAS hits from Nalls 2019 META5 GWAS (https://pdgenetics.shinyapps.io/GWASBrowser/) 
# We will use one as a positive controls for conditional analysis…expect to see little change in effect for these variants
# Include the $ at the end of the rsID to exclude variants that include these rsIDs within their (longer) rsID
grep -e rs7134559$ -e rs10847864$ -e rs11610045$ /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > pos_controls.txt
cat LRRK2_HRC_coding_V2.txt pos_controls.txt > LRRK2_HRC_coding_V3.txt

# Add in the first column ("ID") in CHR:POS format 
cut -f 1,2 LRRK2_HRC_coding_V3.txt | sed 's/\t/:/g' | sed 's/Chr:Start/ID/g' > first_column.txt
paste first_column.txt LRRK2_HRC_coding_V3.txt > LRRK2_HRC_coding_V4.txt

# Copy the IDs back to working directory as LRRK2_coding_VOI.txt (VOI=variants of interest)
# tail -n+2 returns the list without the header
cut -f1 LRRK2_HRC_coding_V4.txt  | tail -n+2 > LRRK2_HRC_coding_V4_IDs.txt
scp LRRK2_HRC_coding_V4_IDs.txt /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt

wc -l LRRK2_HRC_coding_V4_IDs.txt
# 48 variants present…
```

### 3.2 Perform GWAS of chromosome 12 for each cohort

```
# Normal GWAS for IPDGC cohorts
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=24:00:00 IPDGC_CHR12_GWAS.sh NORMAL

# Conditional GWAS for IPDGC cohorts (no rs76904798 + no G2019S)
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=24:00:00 IPDGC_CHR12_GWAS.sh SPECIAL

# Special conditional GWAS for IPDGC cohorts (no N2081D + no G2019S)
sbatch --cpus-per-task=10 --mem=100g --mail-type=ALL --time=24:00:00 IPDGC_CHR12_GWAS.sh CONDI
```

```
# This is IPDGC_CHR12_GWAS.sh

#!/bin/bash

# sh IPDGC_CHR12_GWAS.sh NORMAL
# sh IPDGC_CHR12_GWAS.sh SPECIAL
# sh IPDGC_CHR12_GWAS.sh CONDI

GWAS_TYPE=$1

cd /data/LNG/Julie/Julie_LRRK2_Condi

mkdir ${GWAS_TYPE}_GWAS_CHR12
module load plink/2.0-dev-20191128

cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
	plink2 --bfile HARDCALLS_with_rs10847864 --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--chr 12 \
	--keep /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.$line.txt \
	--pheno-name PHENO_PLINK --covar-variance-standardize \
	--pheno /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.$line.txt \
	--covar /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.$line.txt \
	--covar-name AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_GWAS_CHR12/${GWAS_TYPE}_GWAS_CHR12.$line
done

## EXCEPTIONS => VANCE + MF no age...
plink2 --bfile HARDCALLS_with_rs10847864 --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.VANCE.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.VANCE.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.VANCE.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_GWAS_CHR12/${GWAS_TYPE}_GWAS_CHR12.VANCE

plink2 --bfile HARDCALLS_with_rs10847864 --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.MF.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.MF.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_COVARIATES/LRRK2_condi_covariates_${GWAS_TYPE}.MF.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/${GWAS_TYPE}_GWAS_CHR12/${GWAS_TYPE}_GWAS_CHR12.MF
```

### 3.3 Munging data

```
Files => 
NORMAL_GWAS_CHR12.*.PHENO_PLINK.glm.logistic.hybrid
SPECIAL_GWAS_CHR12.*.PHENOTYPE.glm.logistic.hybrid
CONDI_GWAS_CHR12.*.PHENOTYPE.glm.logistic.hybrid

# This is the header we want
ID REF A1 A1_FREQ beta LOG.OR._SE P

# Loop over to reformat the .hybrid files into .txt files
cd /data/LNG/Julie/Julie_LRRK2_Condi/
cat cohort_file.txt | while read line
do 
  cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/
  Rscript --vanilla /data/LNG/Julie/Julie_LRRK2_Condi/reformat_IPDGC.R NORMAL_GWAS_CHR12.$line.PHENO_PLINK.glm.logistic.hybrid
  cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/
  Rscript --vanilla /data/LNG/Julie/Julie_LRRK2_Condi/reformat_IPDGC.R SPECIAL_GWAS_CHR12.$line.PHENO_PLINK.glm.logistic.hybrid
  cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/
  Rscript --vanilla /data/LNG/Julie/Julie_LRRK2_Condi/reformat_IPDGC.R CONDI_GWAS_CHR12.$line.PHENO_PLINK.glm.logistic.hybrid
done

## Organize the files 
cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/
mkdir prep_files
mv *.log prep_files
mv *.hybrid prep_files

cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/
mkdir prep_files
mv *.log prep_files
mv *.hybrid prep_files

cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/
mkdir prep_files
mv *.log prep_files
mv *.hybrid prep_files
```

```
# This is reformat_IPDGC.R

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla reformat_IPDGC.R $FILENAME
# Rscript --vanilla reformat_IPDGC.R NORMAL_GWAS_CHR12.NEUROX_DBGAP.PHENO_PLINK.glm.logistic.hybrid
FILENAME = args[1]
print(args[1])
print(args[2])
print(FILENAME)
require(dplyr)
require(data.table)
data <- read.table(FILENAME, header = F)
data <- data %>% select(3,4,6,13,19,20,22)
names(data) <- c("ID","REF","A1","A1_FREQ","OR","LOG.OR._SE","P")
data$beta <- log(data$OR)
#reorder the columns and get rid of OR
data <- data[c("ID","REF","A1","A1_FREQ","beta","LOG.OR._SE","P")]
write.table(data, file=paste(sub("PHENO.*", "", FILENAME),"txt", sep=""),quote=F,row.names=F,sep="\t")
```

### 3.4 Pull LRRK2 coding variants from IPDGC CHR12 GWAS

Purpose: After incorporating UKB data, we will make forest plots of LRRK2 coding variants and positive controls.

```
# Make new GWAS results files with only the LRRK2 coding VOI

cd /data/LNG/Julie/Julie_LRRK2_Condi/
cat cohort_file.txt | while read line
do 
  cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12
  # Keep the header to add to VOI files
  head -1 NORMAL_GWAS_CHR12.$line.txt > header.txt
  # Pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt NORMAL_GWAS_CHR12.$line.txt > temp.txt
  # Combine these lines with the header
  cat header.txt temp.txt > NORMAL_GWAS_VOI.$line.txt
  rm header.txt temp.txt
  
  cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12
  # Keep the header to add to VOI files
  head -1 SPECIAL_GWAS_CHR12.$line.txt > header.txt
  # Pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt SPECIAL_GWAS_CHR12.$line.txt > temp.txt
  # Combine these lines with the header
  cat header.txt temp.txt > SPECIAL_GWAS_VOI.$line.txt
  rm header.txt temp.txt
  
  cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12
  # Keep the header to add to VOI files
  head -1 CONDI_GWAS_CHR12.$line.txt > header.txt
  # Pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt CONDI_GWAS_CHR12.$line.txt > temp.txt
  # Combine these lines with the header
  cat header.txt temp.txt > CONDI_GWAS_VOI.$line.txt
  rm header.txt temp.txt
done

## Reorganize the files --> save the VOI files in LRRK2_coding_VOI within each GWAS directory

cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12
mkdir LRRK2_coding_VOI
mv NORMAL_GWAS_VOI* LRRK2_coding_VOI

cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12
mkdir LRRK2_coding_VOI
mv SPECIAL_GWAS_VOI* LRRK2_coding_VOI

cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12
mkdir LRRK2_coding_VOI
mv CONDI_GWAS_VOI* LRRK2_coding_VOI
```

```
##Sanity check: make sure all of the VOI files have 39 lines (including header)
#These are the variants from LRRK2_coding_VOI.txt that are present in our data

cd /data/LNG/Julie/Julie_LRRK2_Condi/
cat cohort_file.txt | while read line
do 
  cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI
  wc -l NORMAL_GWAS_VOI.$line.txt
  cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI
  wc -l SPECIAL_GWAS_VOI.$line.txt
  cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI
  wc -l CONDI_GWAS_VOI.$line.txt
done
```

## 4. Add in UKBiobank
 
This section goes through: 
- Subsetting the UK Biobank data
- Making covariate files
- Performing GWAS on CHR12
- Reformatting data for meta-analysis

### 4.1 Overview of the data

```
# Raw data filtered and unrelated:
/data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.fam 

# PD cases:
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD.txt 

# PD proxy are here:
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_parent_no_PD.txt 

# Samples to NOT use as controls are here:
/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_case_or_PD_parent.txt 

# Select controls (AGE >= 60) from here:
/data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.*
# but make sure you exclude the PD cases and proxies...
# Ideally case-control ratio is 1:10

# Full covariates are here:
/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt
# to include in GWAS are: PC1-5 (to be generated), TOWNSERND, AGE of RECRUIT and SEX

# To create PCs from here:
/data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.*
```

### 4.2 Filter UKB data

```
cd /data/LNG/Julie/Julie_LRRK2_Condi
mkdir UKB_GWAS

#filter for CHR12, Europeans only, etc...
#!/bin/bash
sbatch --cpus-per-task=20 --mem=240g --mail-type=END --time=24:00:00 make_pfile_UKB.sh
module load plink/2.0-dev-20191128
cd /data/CARD/UKBIOBANK/IMPUTED_DATA/
#Extract CHR12 SNPs, Europeans only
plink2 --bgen ukb_imp_chr12_v3.bgen --extract CHR12.SNPS_0_8.txt --geno 0.1 --hwe 1e-6 \
--keep EUROPEAN.txt --make-pgen --mind 0.1 --sample ukb33601_imp_chr1_v3_s487395.sample \
--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/chr12.UKBB.EU.filtered_NEW \
--memory 235000
```

### 4.3 Subset phenotype data for covariate files

Coviariate files to be made: N=6
1. PD vs control 
2. PD vs control (no rs76904798 + no G2019S)
3. PD vs control (no N2081D + no G2019S)
4. PD proxy vs control 
5. PD proxy vs control (no rs76904798 + no G2019S)
6. PD proxy vs control (no N2081D + no G2019S)

Working directory: data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

#### Determine LRRK2 carrier status 

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

module load plink/2.3-alpha
plink2 --bgen /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb_imp_chr12_v3.bgen --snps rs76904798,rs34637584,rs33995883 --make-bed --sample /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb33601_imp_chr1_v3_s487395.sample --out LRRK2_area_snps

module load plink
plink --bfile LRRK2_area_snps --recodeA --out LRRK2_area_snps2

# LRRK2 carrier status stored here: LRRK2_area_snps2.raw

### Checking imputation quality:
# rs76904798 => 5' risk variant
# rs34637584 => G2019S
# rs33995883 => N2081D
# rs10847864 => positive control

grep -e rs76904798 -e rs34637584 -e rs33995883 -e rs10847864 /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb_mfi_chr12_v3.txt

NAME		RS		BP		REF	ALT	MAF		ALT	R2		
12:40614434_C_T	rs76904798	40614434	C	T	0.144409	T	0.996952
rs34637584	rs34637584	40734202	G	A	0.000321086	A	1
rs33995883	rs33995883	40740686	A	G	0.0163887	G	1
rs10847864	rs10847864	123326598	G	T	0.350165	T	1
```

#### Subset phenotype info based on LRRK2 status

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

# Make directory for co-inheritance analysis to copy over some files
mkdir co_inheritance

# Organize some files
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
mkdir LRRK2_status_prep
mv LRRK2_area_snps* LRRK2_status_prep

module load R
R
require(data.table)
require(dplyr)

### Read in all of the files:

# Full covariate file
cov <- fread("/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt",header=T)
# PD cases 
PD <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD.txt",header=T)
# PD proxy cases
Proxy <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_parent_no_PD.txt",header=F)
# Raw data filtered and unrelated -- use to filter the final groups
keep <- fread("/data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.fam",header=F)
# Not controls: don't use PD cases or proxy cases as controls
not_control <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_case_or_PD_parent.txt",header=T)
# Use LRRK2 status to subset covariate files for separate GWAS
LRRK2_status <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/LRRK2_status_prep/LRRK2_area_snps2.raw",header=T)

# Pull the sample IDs for each of the groups
PDshort <-  data.frame(PD$eid)
Proxyshort <- data.frame(Proxy$V1)
not_controlshort <- data.frame(not_control$eid)
keepshort <- data.frame(keep$V1)

# Pull the relevant phenotype info
covshort <- data.frame(cov$FID, cov$AGE_OF_RECRUIT, cov$GENETIC_SEX)

# Rename the sample IDs with "FID" 
names(PDshort) <- c("FID")
names(Proxyshort) <- c("FID")
names(not_controlshort) <- c("FID")
names(keepshort) <- c("FID")

# Rename the cov file header
names(covshort) <- c("FID", "AGE", "SEX")

# Merge the phenotype info with the FIDs for each group to add age and sex info
PD2 = merge(PDshort,covshort,by.x='FID',by.y='FID')
Proxy2 = merge(Proxyshort,covshort,by.x='FID',by.y='FID')
not_control2 = merge(not_controlshort,covshort,by.x='FID',by.y='FID')

# Filter the groups using the filtered and unrelated data (keepshort)
PD3 = merge(PD2,keepshort,by.x='FID',by.y='FID')
Proxy3 = merge(Proxy2,keepshort,by.x='FID',by.y='FID')
not_control3 = merge(not_control2,keepshort,by.x='FID',by.y='FID')

# Get all controls....
keep3 = merge(keepshort,covshort,by.x='FID',by.y='FID')
# Use anti_join to find unmatched records, those that are in keep3 but not not_control3 --> controls 
keep4 = anti_join(keep3,not_control3, by = c('FID'))
# We only want controls that are 60+ yrs old
keep5 = subset(keep4, AGE >= 60)

# Add column on PD status
PD3$STATUS <- "PD"
Proxy3$STATUS <- "PROXY"
keep5$STATUS <- "CONTROL"

### Subset controls randonly... 
# We want 10X as many controls as there are cases
# Case-control = 1530 x 10 = 15300
nrow(PD3) # 1530
# Remember that keep5 is the controls >= 60 yrs
# Set a random seed so you choose the same controls if re-run
set.seed(76418)
PD3_controls <- sample_n(keep5, 15300)
# Use the rest of the controls for the Proxy GWAS
Proxy3_controls <- anti_join(keep5,PD3_controls, by = c('FID'))

# Cat dataframes
PD_FINAL <- rbind(PD3, PD3_controls)
PROXY_FINAL <- rbind(Proxy3, Proxy3_controls)
# Add the IID column
PD_FINAL$IID <- PD_FINAL$FID
PROXY_FINAL$IID <- PROXY_FINAL$FID 
# Rearrange so that IID is the second column and age is last
PD_FINAL <- PD_FINAL[c(1,5,3,4,2)]
PROXY_FINAL <- PROXY_FINAL[c(1,5,3,4,2)]

# Subset based on LRRK2 carrier status
# rs76904798_T => 5' risk variant
# rs34637584_A => G2019S 
# rs33995883_G => N2081D

# Get rid of unneeded columns: only need FID, rs76904798_T, rs34637584_A, rs33995883_G
LRRK2_status$IID <- NULL
LRRK2_status$PAT <- NULL
LRRK2_status$MAT <- NULL
LRRK2_status$SEX <- NULL
LRRK2_status$PHENOTYPE <- NULL
   
# Add LRRK2 status to the final PD and Proxy datasets
PD_FINAL_LRRK2 <- merge(PD_FINAL,LRRK2_status,by.x='FID',by.y='FID')
PROXY_FINAL_LRRK2 <- merge(PROXY_FINAL,LRRK2_status,by.x='FID',by.y='FID')

# Subset the final PD and Proxy datasets based on LRRK2  status
PD_FINAL_norisk_GS <- subset(PD_FINAL_LRRK2, rs76904798_T == 0 & rs34637584_A == 0)
PD_FINAL_noN2081D_GS <- subset(PD_FINAL_LRRK2, rs33995883_G == 0 & rs34637584_A == 0)
PROXY_FINAL_norisk_GS <- subset(PROXY_FINAL_LRRK2, rs76904798_T == 0 & rs34637584_A == 0)
PROXY_FINAL_noN2081D_GS <- subset(PROXY_FINAL_LRRK2, rs33995883_G == 0 & rs34637584_A == 0)

# Add these for co-inheritance analysis later, not for GWAS
PD_FINAL_noGS <- subset(PD_FINAL_LRRK2, rs34637584_A == 0)
PD_FINAL_norisk <- subset(PD_FINAL_LRRK2, rs76904798_T == 0)
PD_FINAL_noND <- subset(PD_FINAL_LRRK2, rs33995883_G == 0)
PROXY_FINAL_noGS <- subset(PROXY_FINAL_LRRK2, rs34637584_A == 0)
PROXY_FINAL_norisk <- subset(PROXY_FINAL_LRRK2, rs76904798_T == 0)
PROXY_FINAL_noND <- subset(PROXY_FINAL_LRRK2, rs33995883_G == 0)

# Save dataframes
# Save PD_FINAL_LRRK2 and PROXY_FINAL_LRRK2 as the normal files so they have the same format as the others
write.table(PD_FINAL_LRRK2, file="UKB_PD_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_norisk_GS, file="UKB_PD_cases_control_over60_noriskGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_noN2081D_GS, file="UKB_PD_cases_control_over60_noNDGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_LRRK2, file="UKB_Proxy_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_norisk_GS, file="UKB_Proxy_cases_control_over60_noriskGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_noN2081D_GS, file="UKB_Proxy_cases_control_over60_noNDGS.txt", quote=FALSE,row.names=F,sep="\t")

# Save these to the co-inheritance directory
write.table(PD_FINAL_noGS, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_PD_cases_control_over60_noGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_norisk, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_PD_cases_control_over60_norisk.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_noND, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_PD_cases_control_over60_noND.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_noGS, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_Proxy_cases_control_over60_noGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_norisk, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_Proxy_cases_control_over60_norisk.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_noND, file="/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/UKB_Proxy_cases_control_over60_noND.txt", quote=FALSE,row.names=F,sep="\t")

q()
n

# Copy some files to the co-inheritance directory
cp UKB_P* /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/
```

### 4.4 Calculate PCs

See here for more information on flashpca: https://github.com/gabraham/flashpca

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
module load flashpca
module load plink

# Test run flashpca on the PD normal group
plink --bfile /data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins \
--maf 0.05 --geno 0.01 --hwe 5e-6 --autosome \
--exclude /data/CARD/GENERAL/exclusion_regions_hg19.txt --make-bed --out FILENAME_2 \
--keep UKB_PD_cases_control_over60.txt
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
flashpca --bfile FILENAME_3 --suffix _UKB_PD_cases_control_over60.txt --numthreads 19
```

#### Run flashpca in a loop

```
# Combine the subset filenames
ls UKB_P* > PC_files.txt

# Make sure PC_files.txt contains these phenptype files: 
# UKB_PD_cases_control_over60_noNDGS.txt
# UKB_PD_cases_control_over60_noriskGS.txt
# UKB_PD_cases_control_over60.txt
# UKB_Proxy_cases_control_over60_noNDGS.txt
# UKB_Proxy_cases_control_over60_noriskGS.txt
# UKB_Proxy_cases_control_over60.txt

sbatch --cpus-per-task=16 --mem=200g --mail-type=ALL --time=24:00:00 flashpca_UKB.sh PC_files.txt
```

```
# This is flashpca_UKB.sh

#!/bin/bash
# sh flashpca_UKB.sh PC_files.txt

PC_files=$1
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

cat $PC_files  | while read line
do 
	plink --bfile /data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins \
	--maf 0.05 --geno 0.01 --hwe 5e-6 --autosome \
	--exclude /data/CARD/GENERAL/exclusion_regions_hg19.txt --make-bed --out FILENAME_2 \
	--keep $line
	plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
	plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
	flashpca --bfile FILENAME_3 --suffix _$line --numthreads 19
done
```

#### Merge the PCs with phenotype info

```
# Merge files in R
module load R
R
require(data.table)
require(dplyr)

# Load the full covariates file
cov <- fread("/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt",header=T)
# Remove the PC columns from cov so that we can merge with the new PCs
cov2 <- cov %>% select(1:8)

# Pull the subset phenptype files from your working directory
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/", pattern="^UKB_P") 
for(file in file.names) {  
  pc <- fread(paste("pcs_",file,sep=""),header=T)
  pc$IID <- NULL
  pheno <- fread(file,header=T)
  # This info will become redundant when merging, keep (PD/CONTROL/PROXY) STATUS and LRRK2 carrier status
  pheno$IID <- pheno$SEX <- pheno$AGE <- NULL
  Mrg = merge(cov2,pheno,by='FID')
  Mrg2 = merge(Mrg,pc,by='FID')
  write.table(Mrg2, file=paste("COV_",file,sep=""), quote=FALSE,row.names=F,sep="\t")
}
q()
n

# Replace PD/CONTROL/PROXY STATUS with numbers
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
cat PC_files.txt | while read line
do 
  sed -i 's/PD/2/g' COV_$line
  sed -i 's/PROXY/2/g' COV_$line
  sed -i 's/CONTROL/1/g' COV_$line
done

# Organize and remove some files
mkdir flashpca_files
mv eigenvalues_UKB_P* flashpca_files
mv eigenvectors_UKB_P* flashpca_files
mv pcs_UKB_P* flashpca_files
mv pve_UKB_P* flashpca_files
rm FILENAME_*
rm pruned_data*
```

### 4.5 Start GWAS on CHR 12

#### Perform GWAS in a loop

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
mkdir GWAS_output
module load plink/2.0-dev-20191128

6 GWAS:
1. COV_UKB_PD_cases_control_over60_noNDGS.txt
2. COV_UKB_PD_cases_control_over60_noriskGS.txt
3. COV_UKB_PD_cases_control_over60.txt
4. COV_UKB_Proxy_cases_control_over60_noNDGS.txt
5. COV_UKB_Proxy_cases_control_over60_noriskGS.txt
6. COV_UKB_Proxy_cases_control_over60.txt

# Test GWAS
# COV_UKB_PD_cases_control_over60_noNDGS.txt
plink2 --pfile chr12.UKBB.EU.filtered_NEW \
--pheno-name STATUS --pheno COV_UKB_PD_cases_control_over60_noNDGS.txt \
--covar COV_UKB_PD_cases_control_over60_noNDGS.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/COV_UKB_PD_cases_control_over60_noNDGS_chr12 --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize

sbatch --cpus-per-task=20 --mem=240g --mail-type=ALL --time=24:00:00 UKB_CHR12_GWAS.sh PC_files.txt
```

```
# This is UKB_CHR12_GWAS.sh

#!/bin/bash
# sh UKB_CHR12_GWAS.sh PC_files.txt

PC_files=$1
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

# Loop for all GWAS chr 12 only
# ${line%%.*} allows you to remove the file extension

cat $PC_files | while read line
do 
	plink2 --pfile chr12.UKBB.EU.filtered_NEW \
	--pheno-name STATUS --pheno COV_$line \
	--covar COV_$line --memory 235000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/COV_${line%%.*}_chr12 --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done
```

#### Check the cases and controls included in each GWAS

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/

# First check the case and control numbers in the covariate files 
R
require(dplyr)
require(data.table)
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/", pattern="^COV_UKB_P") 
for(file in file.names) {  
  cov <- fread(paste("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/",file,sep=""),header=T)
  num_cases <- sum(cov$STATUS==2)
  num_controls <- sum(cov$STATUS==1)
  print(paste(file," has ",num_cases," cases and ",num_controls," controls",sep=""))
}
q()
n

Results: 
[1] "COV_UKB_PD_cases_control_over60_noNDGS.txt has 1466 cases and 14782 controls"
[1] "COV_UKB_PD_cases_control_over60_noriskGS.txt has 1063 cases and 11123 controls"
[1] "COV_UKB_PD_cases_control_over60.txt has 1529 cases and 15279 controls"
[1] "COV_UKB_Proxy_cases_control_over60_noNDGS.txt has 12942 cases and 136307 controls"
[1] "COV_UKB_Proxy_cases_control_over60_noriskGS.txt has 9679 cases and 103040 controls"
[1] "COV_UKB_Proxy_cases_control_over60.txt has 13404 cases and 140655 controls"

# Compare this to the results from the GWAS log files
ls *.log > log_files.txt

cat log_files.txt  | while read line
do 
	echo $line
	grep "1 binary phenotype loaded" $line
done

Results: 
COV_UKB_PD_cases_control_over60_chr12.log
1 binary phenotype loaded (1529 cases, 15279 controls).
COV_UKB_PD_cases_control_over60_noNDGS_chr12.log
1 binary phenotype loaded (1466 cases, 14782 controls).
COV_UKB_PD_cases_control_over60_noriskGS_chr12.log
1 binary phenotype loaded (1063 cases, 11123 controls).
COV_UKB_Proxy_cases_control_over60_chr12.log
1 binary phenotype loaded (13404 cases, 140655 controls).
COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.log
1 binary phenotype loaded (12942 cases, 136307 controls).
COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.log
1 binary phenotype loaded (9679 cases, 103040 controls).
```
| File                                              | COV file |          | GWAS log file |          |
|---------------------------------------------------|----------|----------|---------------|----------|
|                                                   | Cases    | Controls | Cases         | Controls |
| COV_UKB_PD_cases_control_over60_chr12             | 1529     | 15279    | 1529          | 15279    |
| COV_UKB_PD_cases_control_over60_noNDGS_chr12      | 1466     | 14782    | 1466          | 14782    |
| COV_UKB_PD_cases_control_over60_noriskGS_chr12    | 1063     | 11123    | 1063          | 11123    |
| COV_UKB_Proxy_cases_control_over60_chr12          | 13404    | 140655   | 13404         | 140655   |
| COV_UKB_Proxy_cases_control_over60_noNDGS_chr12   | 12942    | 136307   | 12942         | 136307   |
| COV_UKB_Proxy_cases_control_over60_noriskGS_chr12 | 9679     | 103040   | 9679          | 103040   |


### 4.6 Make data ready for meta-analysis

These are the files to process:
- COV_UKB_PD_cases_control_over60_chr12.STATUS.glm.logistic.hybrid
- COV_UKB_PD_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid
- COV_UKB_PD_cases_control_over60_noriskGS_chr12.STATUS.glm.logistic.hybrid
- COV_UKB_Proxy_cases_control_over60_chr12.STATUS.glm.logistic.hybrid
- COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid
- COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.STATUS.glm.logistic.hybrid

#### Reformat plink2 GWAS output
```
# Check the number of variants in these files…
wc -l COV_UKB_PD_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid 
# 1,357,686 variants--in comparison, 439,402 variants in /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr12.UKBB.EU.filtered.pvar

cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output
ls COV_UKB_P*.hybrid | cat > GWAS_files.txt
module load python/3.6

# ${line%%.*} gets rid of the file extension
cat GWAS_files.txt | while read line
do 
	# Filter the results by A1_FREQ (minor allele frequency) >=0.0001 --> to input.txt \
	awk '{ if($13 >= 0.0001) { print }}' $line > input.txt
	
	#reformat the plink2 results for meta-analysis using python
	if [[ $line == *"PD"* ]]; then
		python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
		--outfile toMeta.${line%%.*}.txt --B-or-C B
	elif [[ $line == *"Proxy"* ]]; then
		python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
		--outfile toProxy.${line%%.*}.txt --B-or-C B; fi
done

# The PD files are ready for meta-analysis, the proxy files need more reformatting…
```

```
# This is reformat_plink2_results.py

# -*- coding: utf-8 -*-
"""
# Intro
 - **Project:** GenoML
 - **Author(s):** Mary B. Makarious, Mike A. Nalls, Juan Botia, Hampton Leonard
 - **Date Notebook Started:** 19.09.2019
    - **Quick Description:** Notebook relating to experimental implementation of networkx for community generation.

---
### Quick Description: 
**Problem:** Many software have been used for building networks, within a biological context and in social network and sales analysis. Lets try to do something outside of simple networks.

**Solution:** Lets try to use networkx prioritizing adding Louvain community detection through its widespread application in business https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008. Eventualy we will add decomposition tree based analyses as suggested by https://www.nature.com/articles/s41467-018-03424-4. 

### Motivation/Background:
Data suggests we can do better than very simple coexpression networks also more flexible for SNPs and clinical in nodes?  

### Concerns/Impact on Related Code: 
- Depends on switch to detect discrete or continuous workflows.

### Thoughts for Future Development of Code in this Notebook: 
- Networks could be helpful in partitioning samples and unsupervised learning.

# Imports
"""

import argparse
import sys
import pandas as pd
import numpy as np

"""
# Command args
"""

parser = argparse.ArgumentParser(description='Arguments for training a discrete model')    
parser.add_argument('--infile', type=str, default='plink_glm_file', help='Your results file needing formatting. Default = plink_glm_file.')
parser.add_argument('--outfile', type=str, default='plink_glm_reformatted', help='Your output file. Default = plink_glm_reformatted.')
parser.add_argument('--B-or-C', type=str, default='B', help='B = binary outcome, C = contimuous outcome. Default = B.')
parser.add_argument('--minimal-output', type=str, default='NOPE', help='YEP or NOPE, output only minimal results set to save space. Default = False.')


args = parser.parse_args()

print("")
print("Here is some basic info on the command you are about to run.")
print("Python version info...")
print(sys.version)
print("CLI argument info...")
print("Working with tab-delimited input data", args.infile, "from previous GWAS analyses.")
print("Output file prefix is", args.outfile, "containing your reformatted data in tab delimited format.")
print("Is your GWAS outcome B (binary) or C (continuous)?", args.B_or_C)
print("Only output a minimal set of columns in the results?", args.minimal_output)


print("")

infile = args.infile
outfile = args.outfile
B_or_C = args.B_or_C
minimal_output = args.minimal_output

"""
# Read in your data

This the generic PLINK2 output generated using one of the following commands for binary then continuous outcomes.

module load plink/2.0-dev-20191128

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
plink2 --pfile /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr$chnum.UKBB.EU.filtered \
--pheno-name PHENO_PLINK --pheno covariates/FINAL_cov_UKB_PD_cases_control_over60.txt \
--covar covariates/FINAL_cov_UKB_PD_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out RESULTS/UKB_case_control_chr$chnum --covar-name AGE,SEX,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

or

plink2 --pfile /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr22.UKBB.EU.filtered \
--pheno-name AGE --pheno covariates/FINAL_cov_UKB_Proxy_cases_control_over60.txt \
--covar covariates/FINAL_cov_UKB_Proxy_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out MIKE_continues --covar-name TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize

Just as an additional note, this doesn't scale the contimuous outcome to Z, so be careful with unit standardization across studies with continuous traits.
"""

# For testing
# infile = 'binary.txt'
# outfile = 'test_B_false.tab'
# B_or_C = 'B'
# minimal_output = False

df = pd.read_csv(infile, engine = 'c', sep = '\t')

df.rename(columns = {'#CHROM':'CHROM'}, inplace = True)

print("")
print("Your data looks like this (showing the first few lines of the left-most and right-most columns) ...")
print(df.describe())
print("")
print("Now lets get to the data processing ...")

"""
# Define functions to reformat
"""
def make_log_OR(o):
  b_adjusted = np.log(o)
  return b_adjusted

def reformat_binary_GWAS(x):
    global df_edited
    df_editing = df
    df_editing['markerID'] = df_editing['CHROM'].astype(str) + ":" + df_editing['POS'].astype(str) + ":" + df_editing['REF'].astype(str) + ":" + df_editing['ALT'].astype(str)
    df_editing['effectAllele'] = df_editing['A1']
    df_editing['alternateAllele'] = df_editing['ALT'].where(df_editing['A1'] == df_editing['REF'], df_editing['REF'])
    df_editing['beta']= np.vectorize(make_log_OR)(df_editing['OR'])
    df_editing.rename(columns={'LOG(OR)_SE':'se', 'P':'P', 'A1_FREQ':'effectAlleleFreq', 'OBS_CT':'N', 'ID':'rsID', 'A1_CASE_FREQ':'effectAlleleFreq_cases', 'A1_CTRL_FREQ':'effectAlleleFreq_controls', 'FIRTH?':'firthUsed', 'ERRCODE':'error'}, inplace=True)
    df_edited = df_editing[['markerID','effectAllele','alternateAllele','beta','se','P','effectAlleleFreq','N','rsID','effectAlleleFreq_cases','effectAlleleFreq_controls', 'OR', 'firthUsed','error']]
    del df_editing
    return df_edited

def reformat_continuous_GWAS(x):
    global df_edited
    df_editing = df
    df_editing['markerID'] = df_editing['CHROM'].astype(str) + ":" + df_editing['POS'].astype(str) + ":" + df_editing['REF'].astype(str) + ":" + df_editing['ALT'].astype(str)
    df_editing['effectAllele'] = df_editing['A1']
    df_editing['alternateAllele'] = df_editing['ALT'].where(df_editing['A1'] == df_editing['REF'], df_editing['REF'])
    df_editing.rename(columns={'BETA':'beta', 'SE':'se', 'P':'P', 'A1_FREQ':'effectAlleleFreq', 'OBS_CT':'N', 'ID':'rsID', 'ERRCODE':'error'}, inplace=True)
    df_edited = df_editing[['markerID','effectAllele','alternateAllele','beta','se','P','effectAlleleFreq','N','rsID','error']]
    del df_editing
    return df_edited

if (B_or_C == "B"):
    reformat_binary_GWAS(df)

if (B_or_C == "C"):
    reformat_continuous_GWAS(df)

if (minimal_output == "YEP"):
    df_out = df_edited[['markerID','effectAllele','alternateAllele','beta','se','P','effectAlleleFreq','N']]
    del df_edited

if (minimal_output == "NOPE"):
    df_out = df_edited
    del df_edited

df_out.to_csv(outfile, index=False, sep = '\t')

print("")
print("Here is a little preview of your output (showing the first few lines of the left-most and right-most columns) ...")
print(df_out.describe())
print("")
print("Thanks for using GWAS utilities from GenoML.")

"""
# Quick test
"""
# python reformat_plink2_results.py --infile continuous.txt --outfile test_C_NOPE.tab --B-or-C C --minimal-output NOPE
# python reformat_plink2_results.py --infile continuous.txt --outfile test_C_YEP.tab --B-or-C C --minimal-output YEP
# python reformat_plink2_results.py --infile binary.txt --outfile test_B_NOPE.tab --B-or-C B --minimal-output NOPE
# python reformat_plink2_results.py --infile binary.txt --outfile test_B_YEP.tab --B-or-C B --minimal-output YEP
```

#### Adjust Proxy cases to the same scale as PD cases

```
# Make the toProxy files into .csv files
module load R
R
require(data.table)
data1 <- fread("toProxy.COV_UKB_Proxy_cases_control_over60_chr12.txt",header=T)
data2 <- fread("toProxy.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt",header=T)
data3 <- fread("toProxy.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt",header=T)
write.table(data1, file="toConvert.COV_UKB_Proxy_cases_control_over60_chr12.csv",quote=F,row.names=F,sep=",")
write.table(data2, file="toConvert.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.csv",quote=F,row.names=F,sep=",")
write.table(data3, file="toConvert.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.csv",quote=F,row.names=F,sep=",")
q()
n

# Convert proxies to "normal" files ready for meta-analysis
# This adds the b_adjusted, se_adjusted and p_derived columns
# For meta-analysis in METAL and similar, please use *_adjusted columns. These have been adjusted as per https://www.ncbi.nlm.nih.gov/pubmed/28092683. 
# Taking logistic regression of proxy cases and adjusting to the same scale as actual cases assuming only one parent with disease per proxy case.

python /data/CARD/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
--infile toConvert.COV_UKB_Proxy_cases_control_over60_chr12.csv --beta-proxy beta \
--se-proxy se --p-proxy P --outfile toMeta.COV_UKB_Proxy_cases_control_over60_chr12.csv

python /data/CARD/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
--infile toConvert.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.csv --beta-proxy beta \
--se-proxy se --p-proxy P --outfile toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.csv

python /data/CARD/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
--infile toConvert.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.csv --beta-proxy beta \
--se-proxy se --p-proxy P --outfile toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.csv
```

```
# This is proxy_gwas_gwaxStyle.py

# -*- coding: utf-8 -*-
"""proxy_GWAS_gwaxStyle.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1x-rk7tQvgMzRnX1Rl7WmRRSPtqaTJTUU

# Intro
 - **Project:** proxy_GWAS
 - **Draft author(s):** Mike Nalls and Mary Makarious
 - **Tester(s):** Hampton Leonard and Hirotaka Iwaki
 - **Date Notebook Started:** 16.10.2019
    - **Quick Description:** Notebook to convert proxy GWAS results to GWAS meta-analysis compatible results.

---
### Quick Description: 
**Problem:** Biobanks are a major part of our future work. Although summary stats from case proxies such as 1st degree family history reported individuals are not directly comparable to traditional case-control designs.

**Solution:** Lets adjust these summary stats to be scaled for GWAS meta-analyses. 

### Motivation/Background:
Necessary for GWAS meta-analyses.
See this paper:
https://www.ncbi.nlm.nih.gov/pubmed/28092683

### Concerns/Impact on Related Code: 
- Input summary statistics must be correctly formatted and from a proxy-case analysis.

### Thoughts for Future Development of Code in this Notebook: 
- Would be great to do this on the fly using HAIL [https://hail.is/].

### Testing notes: 
- Tested against http://gwas-browser.nygenome.org/downloads/gwas-browser/gwax-readme.

# Imports
"""

import argparse
import sys
import xgboost
import sklearn
import h5py
import pandas as pd
import numpy as np
import time
import math
import scipy

"""# Insert options for testing in notebook"""

# run_prefix = "test_data.csv"
# beta_name = 'b'
# se_name = 'se'
# p_name = 'p'
# freq_name = 'Freq'
# k_prevalence = 0.2
# scalar = 2
# outfile = "test_data_outfile.csv"

# import os
# from google.colab import drive
# drive.mount('/content/drive/')
# os.chdir("/content/drive/Shared drives/LNG/projects/proxy_GWAS/")
# ! pwd

"""# Command args"""

parser = argparse.ArgumentParser(description='Arguments for training a discrete model')    
parser.add_argument('--infile', type=str, default='GenoML_data', help='Path to your GWAS summary stats data for conversion in comma separated format, should be from logistic regression on family history versus no family history or similar.')
parser.add_argument('--beta-proxy', type=str, default='b', help='Name of beta column.')
parser.add_argument('--se-proxy', type=str, default='se', help='Name of stadnard error column.')
parser.add_argument('--p-proxy', type=str, default='p', help='Name of p vlaue column.')
parser.add_argument('--outfile', type=str, default=1, help='Output CSV path and intended name.')

args = parser.parse_args()

print("")
print("Here is some basic info on the command you are about to run.")
print("Python version info...")
print(sys.version)
print("CLI argument info...")
print("Working with dataset", args.infile, " ... these are youyr GWAS summary stats from the proxy case comparison.")
print("What's your beta coefficient column named?", args.beta_proxy)
print("What's your standard error column named?", args.se_proxy)
print("What's your p value column named?", args.p_proxy)
print("")

run_prefix = args.infile
beta_name = args.beta_proxy
se_name = args.se_proxy
p_name = args.p_proxy
outfile = args.outfile

"""# Read in your data, generally a csv file and summarize"""

df = pd.read_csv(run_prefix, engine = "c")

print("Top few lines fo your file...")
print(df.head())

print("Quick summary of your input data...")
print(df.describe())

print("Quick summary of your data types...")
print(df.dtypes)

"""# Rename columns to b, se, p and freq"""

df.rename(columns={beta_name:'b', se_name:'se', p_name:'p'}, inplace=True)

#df.head(3)

"""# Generate adjusted odds ratio and return adjusted beta"""

print("Rescaling beta.")

def fx_se_scaling(x):
  return x*2

df['b_adjusted'] = np.vectorize(fx_se_scaling)(df['b'])

"""# Generate approximated standard error from adjusted beta and P"""

# Fix floating point errors from small P values in GWAS

print("Now rescaling the SE.")

def fx_b_scaling(x):
  return x*2

df['se_adjusted'] = np.vectorize(fx_b_scaling)(df['se'])

print("Quick sanity check for P derived from adjusted stats. NOTE: \"we keep it real\" regarding floats, so we are capping estimates at 15 decimal places, so there may be some conservative estiamtes but after p < 1E-15 does that really matter?")

def fx_test_P(b, se):
  z_derived = b/se
  p_derived = scipy.special.ndtr(-1 * abs(z_derived)) * 2
  return p_derived

df['p_derived'] = np.vectorize(fx_test_P)(df['b_adjusted'], df['se_adjusted'])

print("For meta-analysis in METAL and similar, please use *_adjusted columns. These have been adjsuted as per https://www.ncbi.nlm.nih.gov/pubmed/28092683. Taking logistic regression of proxy cases and adjusting to the same scale as actual cases assuming only one parent with disease per proxy case.")

"""# Export same dataset with adjusted beta, SE and OR"""

print("Top few lines of your output file...")
print(df.head())

print("Quick summary of your out data...")
print(df.describe())

print("Quick summary of your output data types...")
print(df.dtypes)

df.to_csv(outfile, index=False)
```

```
# Convert the "normalized" proxy .csv files back to .txt files
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output
module load R
R
require(data.table)
data1 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_chr12.csv",header=T)
data2 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.csv",header=T)
data3 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.csv",header=T)
write.table(data1, file="toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt",quote=F,row.names=F,sep="\t")
write.table(data2, file="toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt",quote=F,row.names=F,sep="\t")
write.table(data3, file="toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt",quote=F,row.names=F,sep="\t")
q()
n
```

#### Reformat UKB further and extract LRRK2 coding variants

These are the files to reformat:
- toMeta.COV_UKB_PD_cases_control_over60_chr12.txt
- toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt
- toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt
- toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt
- toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt
- toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt

```
# Copy over the toMeta files for reformatting in a new directory called META
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
mkdir META
cd META
scp /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/toMeta.*.txt .

# We need the variant names of UKB to match those of IPDGC:
# EX of UKB case: 
head -2 /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/META/toMeta.COV_UKB_PD_cases_control_over60_chr12.txt
markerID	effectAllele	alternateAllele	beta	se	P	effectAlleleFreq	N	rsID	effectAlleleFreq_cases	effectAlleleFreq_controls	OR	firthUsed	error
12:87074:T:C	T	C	0.8871296505311677	0.802529	0.268978	0.00041073300000000004	16808	rs564020348	0.0006897289999999999	0.000382813	2.4281.

# Get rid of the :N:N off the end of the markerID for the UKB files
ls | grep UKB > list.txt

cat list.txt
# toMeta.COV_UKB_PD_cases_control_over60_chr12.txt
# toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt
# toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt
# toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt
# toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt
# toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt

cat list.txt  | while read line
do
	sed -i 's/:A//g' $line
	sed -i 's/:T//g' $line
	sed -i 's/:C//g' $line
	sed -i 's/:G//g' $line
done

# These are the headers we want
HEADER of IPDGC:
ID REF A1 A1_FREQ beta LOG.OR._SE P

HEADER of UKB case:
markerID alternateAllele effectAllele effectAlleleFreq beta se P

HEADER of UKB proxy:
markerID alternateAllele effectAllele effectAlleleFreq b_adjusted se_adjusted p_derived
```

```
# Reformat UKB PD files and pull out LRRK2 coding variants
head -1 toMeta.COV_UKB_PD_cases_control_over60_chr12.txt | cut -f 1-7 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > header_PD.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_PD_cases_control_over60_chr12.txt | cut -f 1-7 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > temp
cat header_PD.txt temp > NORMAL_GWAS_VOI.UKBPD.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt | cut -f 1-7 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > temp
cat header_PD.txt temp > SPECIAL_GWAS_VOI.UKBPD.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt | cut -f 1-7 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' > temp
cat header_PD.txt temp > CONDI_GWAS_VOI.UKBPD.txt


# Reformat UKB proxy files and pull out LRRK2 coding variants
# We want to pull the b_adjusted, se_adjusted and p_derived for the proxy cases rather than b, se and p for the PD cases
head -1 toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt | cut -f 1,2,3,7,15,16,17 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > header_proxy.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt | cut -f 1,2,3,7,15,16,17 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > temp
cat header_proxy.txt temp > NORMAL_GWAS_VOI.UKBproxy.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt | cut -f 1,2,3,7,15,16,17 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > temp
cat header_proxy.txt temp > SPECIAL_GWAS_VOI.UKBproxy.txt

grep -f ../../LRRK2_coding_VOI.txt toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt | cut -f 1,2,3,7,15,16,17 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' > temp
cat header_proxy.txt temp > CONDI_GWAS_VOI.UKBproxy.txt


# Copy to folders…
scp NORMAL_GWAS_VOI* /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI
scp SPECIAL_GWAS_VOI* /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI
scp CONDI_GWAS_VOI* /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI
```

#### UK Biobank files to combine with IPDGC for meta-analysis

```
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/NORMAL_GWAS_VOI.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI/SPECIAL_GWAS_VOI.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI/CONDI_GWAS_VOI.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/NORMAL_GWAS_VOI.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI/SPECIAL_GWAS_VOI.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI/CONDI_GWAS_VOI.UKBPD.txt
```

## 5. Making forest plots for LRRK2 coding variants

This section goes through: 
- Pulling the amino acid changes for each variant
- Using metafor to make forest plots for LRRK2 coding variants

#### These are the final headers after reformatting

``` 
HEADER of IPDGC:
ID REF A1 A1_FREQ beta LOG.OR._SE P

HEADER of UKB case:
markerID alternateAllele effectAllele effectAlleleFreq beta se P

HEADER of UKB proxy:
markerID alternateAllele effectAllele effectAlleleFreq b_adjusted se_adjusted p_derived
```

#### Pull the amino acid changes for each variant

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/

module load R
R
require(dplyr)
require(data.table)

# Import some extra info for each of these variants
var_df <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/HRC_LRRK2/LRRK2_HRC_coding_V4.txt",header=T)

# Use the function f to pull out the amino acid change…if not a coding variant, use the rsID
f <- function(row) {
if (row[11] != ".") sub(".*p.", "", row[11]) else row[12]
}
# Note that these names have the one letter amino acid code
AA_short <- c(apply(var_df, 1, f))
id <- var_df$ID

# Split the short amino acid names at the number
split <- do.call(rbind, strsplit(AA_short, "(?<=[A-Z])(?=[0-9])|(?<=[0-9])(?=[A-Z])", perl = TRUE))
first_AA <- split[,1]
middle_num <- split[,2]
last_AA <- split[,3]

# Replace the one letter amino acid code with the three letter code
require(seqinr)
AA_long <- paste(aaa(first_AA),middle_num,aaa(last_AA)) %>% gsub(pattern="NA", replacement="") %>% gsub(pattern=" ",replacement= "")

df <- data.frame(id, AA_short, AA_long)
write.table(df, file="LRRK2_AA_list.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Make individual forest plots for LRRK2 coding variants

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/

make_forest() {
# Create files per variant
gwas_type=$1
cd /data/LNG/Julie/Julie_LRRK2_Condi/${gwas_type}_GWAS_CHR12/LRRK2_coding_VOI
head -1 ${gwas_type}_GWAS_VOI.DUTCH.txt > header.txt
cat /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt | while read line; do
	grep $line ${gwas_type}_GWAS_VOI.*.txt > $line.txt
	cat header.txt $line.txt > header_$line.txt
	sed -e 's/${gwas_type}_GWAS_VOI.//g' header_$line.txt | sed -e 's/.txt:'$line'//g' > header_"$line"v2.txt
done

# Organize the files
mkdir metafor_plots
mv 12:* metafor_plots
mv header* metafor_plots
# Need to copy over the forest plot script so that it can pull the files from its current directory
cp /data/LNG/Julie/Julie_LRRK2_Condi/metafor_LRRK2.R metafor_plots

# Now make the forest plots
cd metafor_plots
cat /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_AA_list.txt | tail -n+2 | while read line; do
	ID=$(echo $line|awk '{print $1}')
	AA=$(echo $line|awk '{print $3}')
	Rscript --vanilla metafor_LRRK2.R $ID $AA ${gwas_type}
done
}

make_forest NORMAL
make_forest CONDI
make_forest SPECIAL

# Copy files
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/*pdf /Users/lakejs/Desktop/NORMAL_forest
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/*pdf /Users/lakejs/Desktop/SPECIAL_forest
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/*pdf /Users/lakejs/Desktop/CONDI_forest

# Make a separate plot for G2019S since it needs different xlim() parameters
cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots
Rscript --vanilla /data/LNG/Julie/Julie_LRRK2_Condi/metafor_LRRK2_G2019S.R 12:40734202 Gly2019Ser 'NORMAL'

scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/12:40734202_NORMAL_final_GS.pdf /Users/lakejs/Desktop
```

```
# This is metafor_LRRK2.R

#!/usr/bin/env Rscript
require(dplyr)
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla metafor_LRRK2.R $FILENAME $FILENAME2 $FILENAME3
# Rscript --vanilla metafor_LRRK2.R 12:40713899 Met1646Thr 'CONDI'
FILENAME = args[1]
FILENAME2 = args[2]
FILENAME3 = args[3]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)
print(FILENAME3)
library(metafor)
data <- read.table(paste("header_",FILENAME,"v2.txt",sep=""), header = T)
## data <- read.table("header_12:40614434v2.txt", header = T)
#pull out the cohort name
labs <- gsub(".*\\.","", data$ID)
yi   <- data$beta
sei  <- data$LOG.OR._SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file = paste(FILENAME,"_",FILENAME3,"_final.pdf",sep=""), width = 8, height = 6)
Pvalue <- formatC(resFe$pval, digits=4)
## pdf(file = "12:40713899_condi_final.pdf", width = 8, height = 6)
forest(resFe, xlim=c(-2,2), main=paste(FILENAME2," P=",Pvalue ,sep=""),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.5, 1, 2, 3)))
dev.off()
```

```
# This is metafor_LRRK2_G2019S.R

#!/usr/bin/env Rscript
require(dplyr)
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla metafor_LRRK2_G2019S.R $FILENAME $FILENAME2 $FILENAME3
# Rscript --vanilla metafor_LRRK2_G2019S.R 12:40734202 Gly2019Ser 'NORMAL'
FILENAME = args[1]
FILENAME2 = args[2]
FILENAME3 = args[3]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)
print(FILENAME3)
library(metafor)
data <- read.table(paste("header_",FILENAME,"v2.txt",sep=""), header = T)
## data <- read.table("header_12:40614434v2.txt", header = T)
# pull out the cohort name
labs <- gsub(".*\\.","", data$ID)
labs <- gsub("NEUROX_DBGAP","NEUROX", labs)
yi   <- data$beta
sei  <- data$LOG.OR._SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file = paste(FILENAME,"_",FILENAME3,"_final_GS.pdf",sep=""), width = 7, height = 6)
## pdf(file = "12:40713899_condi_final.pdf", width = 8, height = 6)
Pvalue <- formatC(resFe$pval, digits=4)
forest(resFe, xlim=c(-4,7),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.25, 0.5, 1, 5, 10, 20, 40)))
mtext(side=3, line = .5, FILENAME2, cex=1.2, font=2)
mtext(side=3, line = -1, paste("P=",Pvalue,sep=""), cex=1, font=2)
dev.off()
```

#### Make combined forest plots with normal, conditional and special conditional GWAS

```
cd /data/LNG/Julie/Julie_LRRK2_Condi
mkdir metafor_combined_plots
cd metafor_combined_plots

# Now make the forest plots
cat /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_AA_list.txt | tail -n+2 | while read line; do
  ID=$(echo $line|awk '{print $1}')
  AA=$(echo $line|awk '{print $3}')
  Rscript --vanilla ../metafor_LRRK2_combined.R $ID $AA
done

# Copy files
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/metafor_combined_plots/*pdf /Users/lakejs/Desktop/Combined_forest
```

```
# This is metafor_LRRK2_combined.R

#!/usr/bin/env Rscript
require(dplyr)
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla metafor_LRRK2_combined.R $FILENAME $FILENAME2
# Rscript --vanilla metafor_LRRK2_combined.R 12:40713899 Met1646Thr
FILENAME = args[1]
FILENAME2 = args[2]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)

library(metafor)
data_normal <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)
data_special <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)
data_condi <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)

labs_normal <- gsub(".*\\.","", data_normal$ID)
#making the plot look better by changing NEUROX_DBGAP to NEUROX
labs_normal <- gsub("NEUROX_DBGAP", "NEUROX", labs_normal)
yi_normal   <- data_normal$beta
sei_normal  <- data_normal$LOG.OR._SE
resFe_normal  <- rma(yi=yi_normal, sei=sei_normal, method="FE")
resRe_normal  <- rma(yi=yi_normal, sei=sei_normal)

labs_special <- gsub(".*\\.","", data_special$ID)
yi_special   <- data_special$beta
sei_special  <- data_special$LOG.OR._SE
resFe_special  <- rma(yi=yi_special, sei=sei_special, method="FE")
resRe_special  <- rma(yi=yi_special, sei=sei_special)

labs_condi <- gsub(".*\\.","", data_condi$ID)
yi_condi   <- data_condi$beta
sei_condi  <- data_condi$LOG.OR._SE
resFe_condi  <- rma(yi=yi_condi, sei=sei_condi, method="FE")
resRe_condi  <- rma(yi=yi_condi, sei=sei_condi)

pdf(file = paste(FILENAME,"_combined.pdf",sep=""), width = 8, height = 7)
Pvalue_normal <- formatC(resFe_normal$pval, digits=4)
Pvalue_special <- formatC(resFe_special$pval, digits=4)
Pvalue_condi <- formatC(resFe_condi$pval, digits=4)

#Make it so that all datasets are included even if NA for the variant
options(na.action = "na.pass")

par(mfrow=c(1,3))

par(mar=c(5,4,1,1))
forest(resFe_normal, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""),
       slab=labs_normal, mlab="Fixed Effects", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1, "Normal", cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_normal,sep=""), cex=1.2, font=2)

par(mar=c(5,0,1,1))
forest(resFe_condi, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=rep("",length(labs_condi)), mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1, "No rs76904798 + No G2019S", cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_condi,sep=""), cex=1.2, font=2)
#adding this for the title
text(0, 18, FILENAME2, cex=1.5, font=2)

par(mar=c(5,0,1,2))
forest(resFe_special, annotate=TRUE, xlim=c(-2.25,3.25), width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=rep("",length(labs_special)), mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5,1, 2, 3)))
text(0, 17.1, "No N2081D + No G2019S", cex=1.2, font=2)
text(0, 16.5, paste("P=",Pvalue_special,sep=""), cex=1.2, font=2)
dev.off()
```

#### Make combined forest plots with normal and special conditional GWAS for rs76904798

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

Rscript --vanilla metafor_combined_rs76904798.R 12:40614434 rs76904798

# Copy the file
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/12:40614434_combined_no_condi.pdf /Users/lakejs/Desktop/
```

```
# This is metafor_combined_rs76904798.R

#!/usr/bin/env Rscript
require(dplyr)
args = commandArgs(trailingOnly=TRUE)
# Start like this
# Rscript --vanilla metafor_combined_rs76904798.R $FILENAME $FILENAME2
# Rscript --vanilla metafor_combined_rs76904798.R 12:40614434 rs76904798
FILENAME = args[1]
FILENAME2 = args[2]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)

library(metafor)
data_normal <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)
data_special <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)

labs_normal <- gsub(".*\\.","", data_normal$ID)
# Making the plot look better by changing NEUROX_DBGAP to NEUROX
labs_normal <- gsub("NEUROX_DBGAP", "NEUROX", labs_normal)
yi_normal   <- data_normal$beta
sei_normal  <- data_normal$LOG.OR._SE
resFe_normal  <- rma(yi=yi_normal, sei=sei_normal, method="FE")
resRe_normal  <- rma(yi=yi_normal, sei=sei_normal)

labs_special <- gsub(".*\\.","", data_special$ID)
yi_special   <- data_special$beta
sei_special  <- data_special$LOG.OR._SE
resFe_special  <- rma(yi=yi_special, sei=sei_special, method="FE")
resRe_special  <- rma(yi=yi_special, sei=sei_special)

pdf(file = paste(FILENAME,"_combined_no_condi.pdf",sep=""), width = 8, height = 7)
Pvalue_normal <- formatC(resFe_normal$pval, digits=4)
Pvalue_special <- formatC(resFe_special$pval, digits=4)

# Make it so that all datasets are included even if NA for the variant
options(na.action = "na.pass")

par(mfrow=c(1,2), oma=c(0,0,2,0))

par(mar=c(5,4,1,1))

forest(resFe_normal, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""),
       slab=labs_normal, mlab="Fixed Effects", col = "red", border = "red", 
       cex=.9, at=log(c(0.5, 1, 2, 3)))
mtext("Normal", line=-1.5, cex=1, font=2)
mtext(paste("P=",Pvalue_normal,sep=""), line=-2.4, cex=1, font=2)

par(mar=c(5,0,1,2))

forest(resFe_special, annotate=TRUE, xlim=c(-2.25,3.25), width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=rep("",length(labs_special)), mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5, 1, 2, 3)))
mtext("No N2081D + No G2019S", line=-1.5, cex=1, font=2)
mtext(paste("P=",Pvalue_special,sep=""), line=-2.4, cex=1, font=2)
mtext(FILENAME2, line=-1, cex=1.4, font=2, outer=TRUE)
dev.off()
```

#### Make combined forest plots with normal and conditional GWAS for N2081D

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

Rscript --vanilla metafor_combined_N2081D.R 12:40740686 Asn2081Asp

# Copy the file
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/12:40740686_combined_no_special.pdf /Users/lakejs/Desktop/
```

```
# This is metafor_combined_N2081D.R

#!/usr/bin/env Rscript
require(dplyr)
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla metafor_combined_N2081D.R $FILENAME $FILENAME2
# Rscript --vanilla metafor_combined_N2081D.R 12:40740686 Asn2081Asp
FILENAME = args[1]
FILENAME2 = args[2]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)

library(metafor)
data_normal <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)
data_condi <- read.table(paste("/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/header_",FILENAME,"v2.txt",sep=""), header = T)

labs_normal <- gsub(".*\\.","", data_normal$ID)
# Making the plot look better by changing NEUROX_DBGAP to NEUROX
labs_normal <- gsub("NEUROX_DBGAP", "NEUROX", labs_normal)
yi_normal   <- data_normal$beta
sei_normal  <- data_normal$LOG.OR._SE
resFe_normal  <- rma(yi=yi_normal, sei=sei_normal, method="FE")
resRe_normal  <- rma(yi=yi_normal, sei=sei_normal)

labs_condi <- gsub(".*\\.","", data_condi$ID)
yi_condi   <- data_condi$beta
sei_condi  <- data_condi$LOG.OR._SE
resFe_condi  <- rma(yi=yi_condi, sei=sei_condi, method="FE")
resRe_condi  <- rma(yi=yi_condi, sei=sei_condi)

pdf(file = paste(FILENAME,"_combined_no_special.pdf",sep=""), width = 8, height = 7)
Pvalue_normal <- formatC(resFe_normal$pval, digits=4)
Pvalue_condi <- formatC(resFe_condi$pval, digits=4)

# Make it so that all datasets are included even if NA for the variant
options(na.action = "na.pass")

par(mfrow=c(1,2), oma=c(0,0,2,0))

par(mar=c(5,4,1,1))

forest(resFe_normal, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""),
       slab=labs_normal, mlab="Fixed Effects", col = "red", border = "red", 
       cex=.9, at=log(c(0.5, 1, 2, 3)))
mtext("Normal", line=-1.5, cex=1, font=2)
mtext(paste("P=",Pvalue_normal,sep=""), line=-2.4, cex=1, font=2)

par(mar=c(5,0,1,2))
forest(resFe_condi, annotate=TRUE, xlim=c(-2.25,3.25),width=3,cex.lab=.8, cex.axis=1,
       atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), 
       slab=rep("",length(labs_condi)), mlab="", col = "red", border = "red", 
       cex=.9, at=log(c(0.5, 1, 2, 3)))
mtext("No rs76904798 + No G2019S", line=-1.5, cex=1, font=2)
mtext(paste("P=",Pvalue_condi,sep=""), line=-2.4, cex=1, font=2)
mtext(FILENAME2, line=-1, cex=1.4, font=2, outer=TRUE)
dev.off()
```

## 6. Check LD co-inheritance of LRRK2 coding variants

This section goes through:
- Making frequency files for IPDGC and UKB data
- Checking co-inheritance of G2019S, rs76904798 and N2081D with all other coding variants
- Determining a frequency cutoff for LRRK2 coding variants

### 6.1 - Make freq files for IPDGC data

#### Make subset files

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load plink

### Create a couple subset files to use

# Recode the genotypes of interest as single allele dosage numbers 
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --snps 12:40734202 --recodeA --out LRRK2_G2019S_only
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --snps 12:40614434 --recodeA --out LRRK2_rs76904798_only
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 --snps 12:40740686 --recodeA --out LRRK2_N2081D_only

module load R
R
data_GS <- read.table("LRRK2_G2019S_only.raw",header=T)
data_rs <- read.table("LRRK2_rs76904798_only.raw",header=T)
data_ND <- read.table("LRRK2_N2081D_only.raw",header=T)

newdata_GS <- subset(data_GS, X12.40734202_A == 0)
newdata_rs <- subset(data_rs, X12.40614434_T == 0) 
newdata_ND <- subset(data_ND, X12.40740686_G == 0) 

dim(newdata_GS) 
# 33757     7
dim(newdata_rs) 
# 33335     7
dim(newdata_ND) 
# 43954     7

# Adding some additional sample info
cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T)

# Drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL

Mrg_GS = merge(newdata_GS,cov,by='FID')
Mrg_rs = merge(newdata_rs,cov,by='FID')
Mrg_ND = merge(newdata_ND,cov,by='FID')

dim(Mrg_GS) 
# 33757    43
dim(Mrg_rs) 
# 33335    43
dim(Mrg_ND) 
# 43954    43

write.table(Mrg_GS, file="LRRK2_G2019S_only_with_COV.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Mrg_rs, file="LRRK2_rs76904798_only_with_COV.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Mrg_ND, file="LRRK2_N2081D_only_with_COV.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

# Excluding all G2019S AND rs76904798 carriers
cp /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_sample_selection.txt .

# Excluding all G2019S AND N2081D carriers
cp /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_special_sample_selection.txt .
```

#### Run association tests in a loop

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load plink

# All data
# Run the normal association tests separately since no file to --keep
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --assoc --out freq
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --logistic --out freq
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --model --out freq

ls LRRK2*.txt > keep_files.txt

cat keep_files.txt
# LRRK2_condi_sample_selection.txt
# LRRK2_condi_special_sample_selection.txt
# LRRK2_G2019S_only_with_COV.txt
# LRRK2_N2081D_only_with_COV.txt
# LRRK2_rs76904798_only_with_COV.txt

# Make a bash array with the output filenames for the conditional tests
# Note that this order corresponds to the order in keep_files.txt
outfile=(
 freq_no_G2019S_rs76904798
 freq_no_G2019S_N2081D
 freq_no_G2019S
 freq_no_N2081D
 freq_no_rs76904798
)

# Run the rest of the association tests 
index=0
cat keep_files.txt | while read line
do
	plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
	--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --assoc --out ${outfile[${index}]} --keep $line
	plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
	--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --logistic --out ${outfile[${index}]} --keep $line
	plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
	--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --model --out ${outfile[${index}]} --keep $line
	index=`expr $index + 1`
done
```

#### Merge files in R

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load R
R
require(dplyr)
require(data.table)

# Make a list of the output files from the association tests
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance", pattern="^freq")

# Pull the prefixes for the output files
file.prefixes <- file.names %>% gsub(pattern="\\..*", replacement="") %>% unique() 

i<-0
dfs = list()
for(file in file.prefixes) { 
 i <- i+1
 data1 <- read.table(paste(file,".assoc",sep=""),header=T)
 data2 <- read.table(paste(file,".assoc.logistic",sep=""),header=T)
 data3 <- read.table(paste(file,".model",sep=""),header=T)
 MM <- data1[,c(2,4,5,6,7)]
 MM2 <- data2[,c(2,6,7,9)]
 MM3 <- subset(data3, TEST=="GENO")
 MM4 <- MM3[,c(2,6,7)]
 MM5 <- merge(MM,MM4,by="SNP")
 MM6 <- merge(MM5,MM2, by="SNP")
 colnames(MM6) <- c("SNP","A1","F_A","F_U","A2","AFF","UNAFF","NMISS","OR","P")
 # Now save the F_A and F_U to the list of dataframes (dfs), for co-inheritance analysis
 filtered <- MM6 %>% select("F_A", "F_U")
 colnames(filtered) <- paste(file, colnames(filtered), sep = "_")
 dfs[[i]]<-data.frame(filtered)
 write.table(MM6, file=paste("LRRK2_coding_variants_", file,".txt",sep=""), quote=FALSE,row.names=F,sep="\t")
}

# Combine the frequency info for each GWAS
SNPs <- MM6$SNP
combined <- bind_cols(SNPs,dfs)

colnames(combined) <- c("SNP","Freq_PD_noNDGS","Freq_Control_noNDGS","Freq_PD_noGSRisk","Freq_Control_noGSRisk","Freq_PD_noGS","Freq_Control_noGS","Freq_PD_noND","Freq_Control_noND","Freq_PD_noRisk","Freq_Control_noRisk","Freq_PD_normal","Freq_Control_normal")

# Reorder based on SNP position, need to use gsub to get rid of the 12: so it sorts properly
SNP_order <- combined$SNP %>% gsub(pattern="12:", replacement="") %>% as.numeric() %>% order()

combined <- combined[SNP_order,]

# Make the header into the first row so it becomes a column with t()
combined <- rbind(colnames(combined), combined)
colnames(combined) <- NULL

# Switch the SNPs to columns
combined <- combined %>% t() %>% data.frame() 

# Set the header to the SNP names
colnames(combined) <- as.character(combined[1,])

# Remove row that's duplicate of header
combined <- combined[-1, ]

write.table(combined, file="IPDGC_freq_tbl.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Make case-control tables for files used in GWAS

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance

# Results including all data (no variant exclusion)
cat LRRK2_coding_variants_freq.txt | cut -f1,3,4,6,7 > case_control_IPDGC_normal.txt

# Results excluding GS and '5 risk variant
cat LRRK2_coding_variants_freq_no_G2019S_rs76904798.txt | cut -f1,3,4,6,7 > case_control_IPDGC_condi.txt

# Results excluding GS and ND variants
cat LRRK2_coding_variants_freq_no_G2019S_N2081D.txt | cut -f1,3,4,6,7 > case_control_IPDGC_special.txt

# Export
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_IPDGC_normal.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_IPDGC_condi.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_IPDGC_special.txt /Users/lakejs/Desktop/
```

### 6.2 Make freq files for UKB data

#### Update the subset files 

```
# Replace PD/CONTROL/PROXY with numbers 
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
ls UKB_P* > UKB_sample_selection.txt

cat UKB_sample_selection.txt
# UKB_PD_cases_control_over60_noGS.txt
# UKB_PD_cases_control_over60_noNDGS.txt
# UKB_PD_cases_control_over60_noND.txt
# UKB_PD_cases_control_over60_noriskGS.txt
# UKB_PD_cases_control_over60_norisk.txt
# UKB_PD_cases_control_over60.txt
# UKB_Proxy_cases_control_over60_noGS.txt
# UKB_Proxy_cases_control_over60_noNDGS.txt
# UKB_Proxy_cases_control_over60_noND.txt
# UKB_Proxy_cases_control_over60_noriskGS.txt
# UKB_Proxy_cases_control_over60_norisk.txt
# UKB_Proxy_cases_control_over60.txt

# Replace PD/CONTROL/PROXY STATUS with numbers
cat UKB_sample_selection.txt | while read line
do 
  sed -i 's/PD/2/g' $line
  sed -i 's/PROXY/2/g' $line
  sed -i 's/CONTROL/1/g' $line
done
```

#### Convert the pfile into bed files using plink2

```
# Make separate bed files since there are separate phenotype files
module load plink/2.3-alpha
cat UKB_sample_selection.txt | while read line
do 
	plink2 --pfile /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/chr12.UKBB.EU.filtered_NEW \
	--keep $line \
	--pheno-name STATUS \
	--pheno $line \
	--make-bed --out ${line%%.*}
done
```

#### Make a new file with the rsIDs of the LRRK2 coding variants

```
# Need the extract file to have the rsIDs to match the UKB data
module load R
R
require(dplyr)
require(data.table)
data <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/HRC_LRRK2/LRRK2_HRC_coding_V4.txt", header=T)

# Substitute an rsID which was missing, found on gnomAD
data[data$avsnp142=="."]$avsnp142 <- "rs768764049"

write.table(data.frame(data$avsnp142,data$ID), file="LRRK2_coding_VOI_rsIDs.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Run association tests in a loop

```
module load plink
cat UKB_sample_selection.txt | while read line
do 
	plink --bfile ${line%%.*} --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --assoc --out ${line%%.*}
	plink --bfile ${line%%.*} --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --logistic --out ${line%%.*}
	plink --bfile ${line%%.*} --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --model --out ${line%%.*}
done
```

#### Merge files in R

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load R
R
require(dplyr)
require(data.table)

# Make a list of the output files from the association tests
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance", pattern="^UKB_P.*cases_control")

# Pull the prefixes for the output files
file.prefixes <- file.names %>% gsub(pattern="\\..*", replacement="") %>% unique() 

# Use to change SNP from the rsID to the CHR:POS format
data <- fread("LRRK2_coding_VOI_rsIDs.txt", header=T)

convert_rsID_to_CHR_POS <- function(SNP) {
 data[data$data.avsnp142 == SNP]$data.ID
}

i<-0
dfs = list()
for(file in file.prefixes) { 
 i <- i+1
 data1 <- read.table(paste(file,".assoc",sep=""),header=T)
 data2 <- read.table(paste(file,".assoc.logistic",sep=""),header=T)
 data3 <- read.table(paste(file,".model",sep=""),header=T)
 MM <- data1[,c(2,4,5,6,7)]
 MM2 <- data2[,c(2,6,7,9)]
 MM3 <- subset(data3, TEST=="GENO")
 MM4 <- MM3[,c(2,6,7)]
 MM5 <- merge(MM,MM4,by="SNP")
 MM6 <- merge(MM5,MM2, by="SNP")
 # Now change the rsID to CHR:POS
 SNPs <-sapply(MM6$SNP, convert_rsID_to_CHR_POS,USE.NAMES=FALSE)
 MM6$SNP <- SNPs
 # Now save the F_A and F_U to the list of dataframes (dfs), for co-inheritance analysis
 filtered <- MM6 %>% select("F_A", "F_U")
 # Use this prefix for the column names to distinguish between groups
 file.prefix <- file %>% gsub(pattern="UKB_|_cases_control_over60|\\.frq.cc", replacement="")
 colnames(filtered) <- paste(file.prefix, colnames(filtered), sep = "_")
 dfs[[i]]<-data.frame(filtered)
 write.table(MM6, file=paste("LRRK2_coding_variants_", file.prefix,".txt",sep=""), quote=FALSE,row.names=F,sep="\t")
}

SNPs <- MM6$SNP
combined <- bind_cols(SNPs,dfs)

# New column names
colnames(combined) <- c("SNP","Freq_PD_noGS","Freq_PDControl_noGS","Freq_PD_noND","Freq_PDControl_noND","Freq_PD_noNDGS","Freq_PDControl_noNDGS","Freq_PD_noRisk","Freq_PDControl_noRisk","Freq_PD_noGSRisk","Freq_PDControl_noGSRisk","Freq_PD_normal","Freq_PDControl_normal","Freq_Proxy_noGS","Freq_ProxyControl_noGS","Freq_Proxy_noND","Freq_ProxyControl_noND","Freq_Proxy_noNDGS","Freq_ProxyControl_noNDGS","Freq_Proxy_noRisk","Freq_ProxyControl_noRisk","Freq_Proxy_noGSRisk","Freq_ProxyControl_noGSRisk","Freq_Proxy_normal","Freq_ProxyControl_normal")

# Reorder based on SNP position, need to use gsub to get rid of the 12: so it sorts properly
SNP_order <- combined$SNP %>% gsub(pattern="12:", replacement="") %>% as.numeric() %>% order()

combined <- combined[SNP_order,]

# Add the header as a row so it's maintained after t()
combined <- rbind(colnames(combined),combined)

# Now delete the header since it's redundant
colnames(combined) <- NULL

# Switch the SNPs to columns
combined <- combined %>% t() %>% data.frame() 

# Set the header to the SNP names
colnames(combined) <- as.character(combined[1,])

# Remove row that's duplicate of header
combined <- combined[-1, ]

write.table(combined, file="UKB_freq_tbl.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Make case-control tables for files used in GWAS

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance

# Results including all data (no variant exclusion)
cat LRRK2_coding_variants_PD.txt | cut -f1,3,4,6,7 > case_control_UKB_PD_normal.txt
cat LRRK2_coding_variants_Proxy.txt | cut -f1,3,4,6,7 > case_control_UKB_Proxy_normal.txt

# Results excluding GS and '5 risk variant
cat LRRK2_coding_variants_PD_noriskGS.txt | cut -f1,3,4,6,7 > case_control_UKB_PD_condi.txt
cat LRRK2_coding_variants_Proxy_noriskGS.txt | cut -f1,3,4,6,7 > case_control_UKB_Proxy_condi.txt

# Results excluding GS and ND variants
cat LRRK2_coding_variants_PD_noNDGS.txt | cut -f1,3,4,6,7 > case_control_UKB_PD_special.txt
cat LRRK2_coding_variants_Proxy_noNDGS.txt | cut -f1,3,4,6,7 > case_control_UKB_Proxy_special.txt

# Export
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_PD_normal.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_Proxy_normal.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_PD_condi.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_Proxy_condi.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_PD_special.txt /Users/lakejs/Desktop/
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/case_control_UKB_Proxy_special.txt /Users/lakejs/Desktop/
```

### 6.3 Determine a freq or count cutoff for the LRRK2 coding variants

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

# Based on the forest plots, should keep these variants
grep -e id -e rs76904798 -e Asn2081Asp -e Leu119Pro -e Asn551Lys -e Ile723Val -e Arg1398His -e Arg1514Gln -e Pro1542Ser -e Ser1647Thr -e Met1646Thr -e Met2397Thr -e rs10847864 LRRK2_AA_list.txt > keep_LRRK2_variants.txt

cat keep_LRRK2_variants.txt
# id	AA_short	AA_long
# 12:40629436	L119P	Leu119Pro
# 12:40657700	N551K	Asn551Lys
# 12:40671989	I723V	Ile723Val
# 12:40702911	R1398H	Arg1398His
# 12:40707778	R1514Q	Arg1514Gln
# 12:40707861	P1542S	Pro1542Ser
# 12:40713899	M1646T	Met1646Thr
# 12:40713901	S1647T	Ser1647Thr
# 12:40740686	N2081D	Asn2081Asp
# 12:40758652	M2397T	Met2397Thr
# 12:40614434	rs76904798	rs76904798
# 12:123326598	rs10847864	rs10847864
```

#### Calculate MAF and minor allele counts for the LRRK2 variants 

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load plink

# IPDGC
plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --freq --out IPDGC_freq

plink --bfile /data/LNG/Julie/Julie_LRRK2_Condi/HARDCALLS_with_rs10847864 \
--extract /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt --freq counts --out IPDGC_freq

# UKB PD
plink --bfile UKB_PD_cases_control_over60 \
 --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --freq --out UKB_PD_freq

plink --bfile UKB_PD_cases_control_over60 \
 --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --freq counts --out UKB_PD_freq

# UKB Proxy
plink --bfile UKB_Proxy_cases_control_over60 \
 --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --freq --out UKB_Proxy_freq

plink --bfile UKB_Proxy_cases_control_over60 \
 --extract <(cut -f1 LRRK2_coding_VOI_rsIDs.txt) --freq counts --out UKB_Proxy_freq
```

#### Merge the frequency and counts info in R

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance
module load R
R
require(dplyr)
require(data.table)

# Import the files
IPDGC_freq <- fread("IPDGC_freq.frq",header=T)
IPDGC_counts <- fread("IPDGC_freq.frq.counts",header=T)
PD_freq <- fread("UKB_PD_freq.frq",header=T)
PD_counts <- fread("UKB_PD_freq.frq.counts",header=T)
Proxy_freq <- fread("UKB_Proxy_freq.frq",header=T)
Proxy_counts <- fread("UKB_Proxy_freq.frq.counts",header=T)
data <- fread("LRRK2_coding_VOI_rsIDs.txt",header=T)
keep <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/keep_LRRK2_variants.txt",header=T)

# Merge the counts and freqs dataframes
PD <- merge(PD_freq,PD_counts,by="SNP")
Proxy <- merge(Proxy_freq,Proxy_counts,by="SNP")
IPDGC <- merge(IPDGC_freq,IPDGC_counts,by="SNP")

# Rename the SNP column to rsID for UKB
PD <- PD %>% rename(rsID = SNP)
Proxy <- Proxy %>% rename(rsID = SNP)

# Rename the ID column to "SNP"
data <- data %>% rename(SNP = data.ID)

# Add the SNP in CHR:POS format to the UKB data and call this column SNP to match IPDGC
PD  <- merge(PD, data, by.x="rsID",by.y="data.avsnp142")
Proxy  <- merge(Proxy, data, by.x="rsID",by.y="data.avsnp142")

# Determine the minor allele count for each of the variants
# This is either C1 or C2 depending on whether A1.y (from _counts) is the minor allele (it is if it matches A1.x from _freq)
# Returns the count of the minor allele
f <- function(row) {if (row["A1.x"] == row["A1.y"]) as.numeric(row["C1"]) else as.numeric(row["C2"])}

PD$Minor_allele_count_UKB_PD <- c(apply(PD, 1, f))
Proxy$Minor_allele_count_UKB_Proxy <- c(apply(Proxy, 1, f))
IPDGC$Minor_allele_count_IPDGC <- c(apply(IPDGC, 1, f))

# Rename the columns 
IPDGC <- IPDGC %>% rename(MAF_IPDGC = MAF) %>% rename(Allele_count_IPDGC = NCHROBS)
PD <- PD %>% rename(MAF_UKB_PD = MAF) %>% rename(Allele_count_UKB_PD = NCHROBS)
Proxy <- Proxy %>% rename(MAF_UKB_Proxy = MAF) %>% rename(Allele_count_UKB_Proxy = NCHROBS)

# Merge all of the dataframes
library(tidyverse)
data2 <- list(IPDGC, PD, Proxy) %>% reduce(inner_join, by = "SNP")

data2 <- data2 %>% select("SNP","MAF_IPDGC","MAF_UKB_PD","MAF_UKB_Proxy","Minor_allele_count_IPDGC","Minor_allele_count_UKB_PD","Minor_allele_count_UKB_Proxy","Allele_count_IPDGC", "Allele_count_UKB_PD","Allele_count_UKB_Proxy")

# Add a column for total allele count and total minor allele count in both datasets combined
data2$Total_allele_count<- data2$Allele_count_IPDGC + data2$Allele_count_UKB_PD + data2$Allele_count_UKB_Proxy
data2$Total_minor_allele_count <- data2$Minor_allele_count_IPDGC + data2$Minor_allele_count_UKB_PD + data2$Minor_allele_count_UKB_Proxy
```

#### Determine a cutoff for which variants to include in final tables

```
## Investingating MAF cutoffs: 
MAF_df <- data2 %>% select("SNP","MAF_IPDGC","MAF_UKB_PD","MAF_UKB_Proxy")

# First see if the variants I want to keep are included if MAF > 0.01 for all three datasets
MAF_0.01 <- MAF_df %>% filter(MAF_IPDGC > 0.01) %>% filter(MAF_UKB_PD > 0.01) %>% filter(MAF_UKB_Proxy > 0.01) 
keep$id %in% MAF_0.01$SNP
# [1] FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

# See what the minimum MAF is for the variants I want to keep
keep_MAF <- MAF_df %>% filter(SNP %in% keep$id)

min(keep_MAF[,2],keep_MAF[,3],keep_MAF[,4])
# [1] 0.002434

# Now try MAF > 0.001 and see what variants are included
MAF_0.001 <- MAF_df %>% filter(MAF_IPDGC > 0.001) %>% filter(MAF_UKB_PD > 0.001) %>% filter(MAF_UKB_Proxy > 0.001) 
keep$id %in% MAF_0.001$SNP
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
 
# Let's see what amino acid changes these SNPs correspond to 
AA <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_AA_list.txt",header=T)

LRRK2_AA_filtered <- AA %>% filter(id %in% MAF_0.001$SNP)

LRRK2_AA_filtered
#              id   AA_short    AA_long
# 1:  12:40629436      L119P  Leu119Pro
# 2:  12:40657700      N551K  Asn551Lys
# 3:  12:40671989      I723V  Ile723Val
# 4:  12:40702911     R1398H Arg1398His
# 5:  12:40707778     R1514Q Arg1514Gln
# 6:  12:40707861     P1542S Pro1542Ser
# 7:  12:40713899     M1646T Met1646Thr
# 8:  12:40713901     S1647T Ser1647Thr
# 9:  12:40740686     N2081D Asn2081Asp
# 10:  12:40758652     M2397T Met2397Thr
# 11:  12:40614434 rs76904798 rs76904798
# 12:  12:46419086  rs7134559  rs7134559
# 13: 12:123326598 rs10847864 rs10847864
# 14: 12:133063768 rs11610045 rs11610045

# The cutoff MAF > 0.001 seems good...

# Add back G2019S and get rid of the two positive controls we won't use (rs7134559 and rs11610045)
GS <- AA[AA$AA_short == "G2019S"]
LRRK2_AA_final <- LRRK2_AA_filtered %>% subset(AA_short != "rs7134559") %>% subset(AA_short != "rs11610045") %>% rbind(GS)

write.table(LRRK2_AA_final, file="LRRK2_AA_final.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

## 7. Make final tables and figures

This section goes through:
- Preparing tables for manuscript
- Preparing figures for manuscript

#### Make a final table with variant info

```
cd /data/LNG/Julie/Julie_LRRK2_Condi
R 
require(dplyr)
require(data.table)

data <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/HRC_LRRK2/LRRK2_HRC_coding_V4.txt",header=T)
data2 <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/LRRK2_AA_final.txt",header=T)
data3 <- merge(data,data2,by.x="ID",by.y="id")
data4 <- data3 %>% select("avsnp142","Chr","Start","Ref","Alt","AA_long")
colnames(data4) <- c("rsID","Chr","Position","Ref","Alt","Amino Acid Change")
data4 <- data4[order(data4$Position)]

# Get rid of the rsIDs in the Amino Acid Change column 
f <- function(row) {
if (startsWith(row["Amino Acid Change"], "rs")) sub(".*", "-", row["Amino Acid Change"]) else row["Amino Acid Change"]
}

AA_change_final <- c(apply(data4, 1, f))
data4$"Amino Acid Change" <- AA_change_final

write.table(data4, file="LRRK2_variant_info.txt", quote=FALSE,row.names=F,sep="\t")

q()
n

# Copy the file
scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_variant_info.txt /Users/lakejs/Desktop/
```

#### Organize the final forest plots into a new directory

```
cd /data/LNG/Julie/Julie_LRRK2_Condi
mkdir final_plots
cd final_plots

module load R 
R
require(data.table)
require(dplyr)

AA_final <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/co_inheritance/LRRK2_AA_final.txt",header=T)

# First make a list of the IDs of the combined plots that we want to keep
dont_keep = c("G2019S","N2081D","rs76904798")
`%notin%` <- Negate(`%in%`)
IDs <- AA_final %>% subset(AA_short %notin% dont_keep) %>% select(id) %>% unlist(use.names=FALSE)

# Make a list of the filenames of all of the combined plots
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/metafor_combined_plots", pattern="^12:")

# Add the full path to these filenames
file.names <- paste("/data/LNG/Julie/Julie_LRRK2_Condi/metafor_combined_plots/", file.names, sep="")

# Now filter file.names for the plots we want to keep
require(sjmisc)
bool <- sapply(file.names, str_contains, pattern=IDs,logic="or",USE.NAMES=FALSE)
combined_files <- file.names[bool]

# Make a list of all of the pdf files to put into the final_plots directory
# Add the three plots that don't have normal, conditional and special GWAS results 
other_files = c(
"/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/LRRK2_coding_VOI/metafor_plots/12:40734202_NORMAL_final_GS.pdf",
"/data/LNG/Julie/Julie_LRRK2_Condi/12:40614434_combined_no_condi.pdf",
"/data/LNG/Julie/Julie_LRRK2_Condi/12:40740686_combined_no_special.pdf"
)

all_files <- append(combined_files, other_files)
write.table(all_files, file="forest_plot_filenames.txt", quote=FALSE,row.names=F,sep="\t",col.names = F)

q()
n

# Copy the files to the current directory final_plots
cat forest_plot_filenames.txt | xargs -i scp {} .

scp lakejs@biowulf.nih.gov://data/LNG/Julie/Julie_LRRK2_Condi/final_plots/*.pdf /Users/lakejs/Desktop/final_plots
```

## Done....

![myImage](https://media.giphy.com/media/XRB1uf2F9bGOA/giphy.gif)


