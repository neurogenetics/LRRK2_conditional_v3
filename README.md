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
3) Perform GWAS excluding risk and G2019S variant on a cohort level basis and meta-analyze...

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

### [3. Perform cohort-level GWAS on IPDGC data excluding 5' risk (rs76904798) and G2019S variants](#3-Perform-cohort-level-GWAS-on-IPDGC-data-excluding-risk-(rs76904798)-and-G2019S-variants)
This section goes through: 
- Extracting all LRRK2 coding variants
- Performing GWAS of chromosome 12 for each cohort
- Prepping before meta-analysis

### [4. Adding in UKBiobank](#4-Adding-in-UKBiobank)
This section goes through: 
- Adding the UK biobank data to point 2

### [5. Combining all data together](#5-Combining-all-data-together)

### [6. Check LD co-inheritance of LRRK2 coding variants](#6-Check-LD-co-inheritance-of-LRRK2-coding-variants)
This section goes through:
- Checking if LRRK2 variants are typically co-inherited or not?
- Assessing frequency of these variants in full data-set
- LRRK2 G2019S with all other coding variants
- LRRK2 rs76904798 with all other coding variants
- Preparing files for Tables for manuscript

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
Variants of interest:
LRRK2 G2019S => hg19 12:40734202:G:A
rs76904798 => hg19 12:40614434:C:T

cd /data/LNG/Julie/Julie_LRRK2_Condi 
module load plink

# simple test
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202,12:40614434 --assoc --out test
# Among remaining phenotypes, 21478 are cases and 24388 are controls.
 CHR           SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR 
  12   12:40614434   40614434    T   0.1564   0.1406    C        45.29    1.702e-11        1.133 
  12   12:40734202   40734202    A 0.007297 0.0005648    G        203.6    3.346e-46        13.01 
# confirming associations, so thats good :)
```

### 1.2 - Checking the imputation quality of the IPDGC data

```
cd /data/LNG/CORNELIS_TEMP/PD_AAO/IMPUTATION_QUALITY/
grep 12:40614434 info_all.12 > LRRK2_1.txt
grep 12:40734202 info_all.12 > LRRK2_2.txt
cat header LRRK2_1.txt LRRK2_2.txt > overview_for_LRRK2_conditional_GWAS.txt

# Copy file home
scp lakejs@biowulf.nih.gov://data/LNG/CORNELIS_TEMP/PD_AAO/IMPUTATION_QUALITY/overview_for_LRRK2_conditional_GWAS.txt /Users/lakejs/Desktop/

# summary
12:40614434 -> very well imputed, present in almost all data
12:40734202 -> also not bad...
```
#### IPDGC imputation quality table:

| SNP	| 12:40614434 | How | 12:40734202 | How |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| DUTCH_Rsq  | 0.99903  | Imputed  | 0.87355  | Imputed |
| FINLAND_Rsq  | 0.99972  | Imputed  | **0.00109**  | Imputed |
| GERMANY_Rsq  | 0.9961  | Imputed  | 0.85943  | Imputed |
| HBS_Rsq  | 0.99999  | Genotyped  | 0.99982  | Genotyped |
| MCGILL_Rsq  | 1  | Genotyped  | 0.91193  | Imputed |
| MF_Rsq  | 0.99897  | Imputed  | **0.67548**  | Imputed |
| NEUROX_DBGAP_Rsq  | 0.99684  | Imputed  | 0.9629  | Imputed |
| NIA_Rsq  | 0.99998  | Genotyped  | 0.98795  | Genotyped |
| OSLO_Rsq  | 0.99775  | Imputed  | 0.93746  | Imputed |
| PDBP_Rsq  | 0.99942  | Imputed  | **0.75988**  | Imputed |
| PPMI_Rsq  | 0.99999  | Genotyped  | 0.97688  | Genotyped |
| SHULMAN_Rsq  | 0.99482  | Imputed  | **0.00335**  | Imputed |
| SPAIN_Rsq  | 0.99557  | Genotyped  | 0.99997  | Genotyped |
| SPAIN2_Rsq  | 0.988  | Imputed  | 0.99957  | Genotyped |
| UK_GWAS_Rsq  | 1  | Genotyped  | 0.99984  | Genotyped |
| VANCE_Rsq  | 1  | Genotyped  | 0.99384  | Imputed |

Also see: LNG G-suite --> users/leonardhl/LRRK2_conditional/Imputation_quality.xlsx

### 1.3 - Assessing frequency of LRRK2 G2019S and rs76904798 in the data
```
# check allelic distribution
cd /data/LNG/Julie/Julie_LRRK2_Condi
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202,12:40614434 --model --out allelic_dist

# allelic distribution:
 CHR           SNP   A1   A2     TEST            AFF          UNAFF        CHISQ   DF            P
  12   12:40614434    T    C     GENO 558/5602/15318 486/5885/18017        46.02    2    1.014e-10
  12   12:40734202    A    G     GENO    1/236/16072     0/20/17685           NA   NA           NA

# so in total to use likely ~25,000 samples because we want homozygous reference for both variants
```
Also see: LNG G-suite --> users/leonardhl/LRRK2_conditional/allelic_dist.xlsx

### 1.4 - Overview of the full data and selection of data

```
##recode the genotypes of interest as single allele dosage numbers 

#use for the no rs76904798 + no G2019S dataset
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202,12:40614434 --recodeA --out LRRK2_condi_variant_selection

#use for the no N2081D + no G2019S dataset
#call this the "special" conditional GWAS: see if rs76904798 signal remains
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202,12:40740686 --recodeA --out LRRK2_condi_special_variant_selection
```

```
##subset the data to include only homozygous reference carriers of both variants

module load R
R
data <- read.table("LRRK2_condi_variant_selection.raw",header=T)
data2 <- read.table("LRRK2_condi_special_variant_selection.raw",header=T)

#X12.40614434_T is rs76904798 and X12.40734202_A is G2019S
newdata <- subset(data, X12.40614434_T == 0 & X12.40734202_A == 0) 
dim(newdata) # 24532     8

#X12:40740686_G is N2081D and X12.40734202_A is G2019S
newdata2 <- subset(data2, X12.40740686_G == 0 & X12.40734202_A == 0)
dim(newdata2) # 32227     8

#add some additional sample info
cov <- read.table("/data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/IPDGC_all_samples_covariates.txt",header=T)
# drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL

MM = merge(newdata,cov,by='FID')
dim(MM) # 24532    44
MM2 = merge(newdata2,cov,by='FID')
dim(MM2) # 32227    44

#datasets with good data (the ones with homo-ref carriers):
# [1] "DUTCH"        "GERMANY"      "HBS"          "MCGILL"       "MF"          
# [6] "NEUROX_DBGAP" "NIA"          "PDBP"         "PPMI"         "SHULMAN"     
# [11] "SPAIN3"       "SPAIN4"       "VANCE"       

#create a file for selecting individuals who are homo-ref carriers
write.table(MM, file="LRRK2_condi_sample_selection.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM2, file="LRRK2_condi_special_sample_selection.txt", quote=FALSE,row.names=F,sep="\t")

#display the case-control distribution for each dataset
require(dplyr)
MM_grouped <- MM %>% group_by(DATASET) %>% summarise(Case = sum(PHENOTYPE == 2), Control = sum(PHENOTYPE == 1), TOTAL = n()) %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "SUM"))
write.table(MM_grouped, file="LRRK2_condi_sample_selection_grouped.txt", quote=FALSE,row.names=F,sep="\t")

MM2_grouped <- MM2 %>% group_by(DATASET) %>% summarise(Case = sum(PHENOTYPE == 2), Control = sum(PHENOTYPE == 1), TOTAL = n()) %>% bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "SUM"))
write.table(MM2_grouped, file="LRRK2_condi_special_sample_selection_grouped.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

#copy files
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

Create cohort loop over file:
/data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt

Create sample inclusion file:
Using LRRK2_condi_sample_selection.txt from above...

Also note that most data was already preprocessed according to this https://github.com/neurogenetics/GWAS-pipeline 
and has been used in previous GWAS such as https://pubmed.ncbi.nlm.nih.gov/31701892/ and https://pubmed.ncbi.nlm.nih.gov/30957308/



### 2.1 Create new PC's for each cohort

First loop over all cleaned unimputed data to create fresh PC's

#### PC files without G2019S and rs76904798

```
### loop for making PCA
cd /data/LNG/Julie/Julie_LRRK2_Condi
cat cohort_file.txt | while read line
do 
	#first copy over the pre_impute_vcf_files to new directories
	mkdir $line
	cp /data/LNG/CORNELIS_TEMP/PD_AAO/pre_impute_vcf_files/$line/{$line.bed,$line.bim,$line.fam} /data/LNG/Julie/Julie_LRRK2_Condi/$line
	
	cd $line
	plink --bfile $line --keep /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_sample_selection.txt --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_CONDI
	
	#send the .eigenvec files back to the working directory to combine into a new combined PC file
	scp $line.LRRK2_condi_PCA.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	cd ..
done

```


#### PC files without G2019S and N2081D

```
### loop for making PCA
cat cohort_file.txt | while read line
do 
	cd $line
	plink --bfile $line --keep /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_condi_special_sample_selection.txt --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_SPECIAL
	
	#send the .eigenvec files back to the working directory to combine into a new combined PC file
	scp $line.LRRK2_condi_PCA_SPECIAL.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	cd ..
done

```

#### PC files with G2019S, N2081D and rs76904798 <= meaning normal files...

```
## loop for making PCA
#all samples!
cat cohort_file.txt  | while read line
do 
	cd $line
	#the only difference here is that we don't use --keep to filter out people with the variants
	plink --bfile $line --maf 0.01 --geno 0.15 --hwe 1E-6 --make-bed --out $line.filter
	plink --bfile $line.filter --indep-pairwise 50 5 0.5 --out prune
	plink --bfile $line.filter --extract prune.prune.in --make-bed --out prune 
	plink --bfile prune --pca --out $line.LRRK2_condi_PCA_NORMAL
	
	#send the .eigenvec files back to the working directory to combine into a new combined PC file
	scp $line.LRRK2_condi_PCA_NORMAL.eigenvec /data/LNG/Julie/Julie_LRRK2_Condi
	cd ..
done

```

### 2.2 Create covariate files

#### Merge new PC's in R with other phenotype data

```
#combine the eigenvec files for each GWAS
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

#we want to keep the phenotype information but get rid of the old PCs before merging
NORMAL_cov2 <- NORMAL_cov[,c(1:10)]
SPECIAL_cov2 <- SPECIAL_cov[,c(1:14)]
CONDI_cov2 <- condi_cov[,c(1:14)]

#now add the new PCs
NORMAL_Mrg <- merge(NORMAL_cov2,NORMAL_PC,by.x="FID",by.y="V1") 
SPECIAL_Mrg <- merge(SPECIAL_cov2,SPECIAL_PC,by.x="FID",by.y="V1") 
CONDI_Mrg <- merge(CONDI_cov2,CONDI_PC,by.x="FID",by.y="V1") 

#get rid of this column since it is a repeat of "FID"
NORMAL_Mrg$V2 <- NULL 
SPECIAL_Mrg$V2 <- NULL 
CONDI_Mrg$V2 <- NULL 

#only keep the first 10 PCs
NORMAL_Mrg2 <- NORMAL_Mrg[,c(1:20)]
SPECIAL_Mrg2 <- SPECIAL_Mrg[,c(1:24)]
CONDI_Mrg2 <- CONDI_Mrg[,c(1:24)]

#change the name of the first 10 PCs

#do this for the normal files
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

#do this for the special conditional files
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

#do this for the conditional files
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

#export the final coviariate files
write.table(NORMAL_Mrg2, file="LRRK2_condi_covariates_NORMAL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(SPECIAL_Mrg2, file="LRRK2_condi_covariates_SPECIAL.txt", quote=FALSE,row.names=F,sep="\t")
write.table(CONDI_Mrg2, file="LRRK2_condi_covariates.txt", quote=FALSE,row.names=F,sep="\t")

q()
n
```

#### Subset each cohort due to potentially overlapping sample names

```
cd /data/LNG/Julie/Julie_LRRK2_Condi

#keep the header so that you can add it to each separate file
head -1 LRRK2_condi_covariates_NORMAL.txt > header1.txt
head -1 LRRK2_condi_covariates_SPECIAL.txt > header2.txt
head -1 LRRK2_condi_covariates.txt > header3.txt

cat cohort_file.txt  | while read line
do
	#pull out all of the lines that contain the cohort name 
	grep $line LRRK2_condi_covariates_NORMAL.txt > temp1
	grep $line LRRK2_condi_covariates_SPECIAL.txt > temp2
	grep $line LRRK2_condi_covariates.txt > temp3
	#combine these lines with the header
	cat header1.txt temp1 > LRRK2_condi_covariates_NORMAL.$line.txt
	cat header2.txt temp2 > LRRK2_condi_covariates_SPECIAL.$line.txt
	cat header3.txt temp3 > LRRK2_condi_covariates.$line.txt
done


# fix MF data...
#the MF covariates file also pulled out other cohorts since "MF" is contained within other cohorts' IDs (e.g. HBS_PD_INVDG562MF3)

#want to return all non matching lines -- lines without HBS, PDBP, SPAIN4
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_NORMAL.MF.txt > temp
mv temp LRRK2_condi_covariates_NORMAL.MF.txt
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates_SPECIAL.MF.txt > temp
mv temp LRRK2_condi_covariates_SPECIAL.MF.txt
grep -v -e HBS -e PDBP -e SPAIN4 LRRK2_condi_covariates.MF.txt > temp
mv temp LRRK2_condi_covariates.MF.txt

# fixed…

```

#### Reorganize the files

```
#move the NORMAL covariate files into a new directory
mkdir NORMAL_COVARIATES
mv *NORMAL.eigenvec NORMAL_COVARIATES
mv LRRK2_condi_covariates_NORMAL* NORMAL_COVARIATES

#move the SPECIAL covariate files into a new directory
mkdir SPECIAL_COVARIATES
mv *SPECIAL.eigenvec SPECIAL_COVARIATES
mv LRRK2_condi_covariates_SPECIAL* SPECIAL_COVARIATES

#move the conditional covariate files into a new directory
mkdir CONDI_COVARIATES
mv *.eigenvec CONDI_COVARIATES
mv LRRK2_condi_covariates* CONDI_COVARIATES

```


## 3. Perform cohort-level GWAS on IPDGC data excluding 5' risk (rs76904798) and G2019S variants

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

#pull out the variants from CHR12, then filter by LRRK2, exonic, nonsynonymous --> LRRK2 coding variants
awk '$1 == "12"' /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt | grep LRRK2 | grep exonic | sed 's/ /_/g' | grep nonsynonymous > temp
cat header.txt temp > LRRK2_HRC_coding.txt

# manually add in rs76904798 (non-coding)
grep rs76904798 /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > risk_variant.txt
cat LRRK2_HRC_coding.txt risk_variant.txt > LRRK2_HRC_coding_V2.txt

#manually add in non-LRRK2 CHR12 GWAS hits from Nalls 2019 META5 GWAS (https://pdgenetics.shinyapps.io/GWASBrowser/) as positive controls for conditional analysis…expect to see little change in effect for these variants 
grep -e rs7134559$ -e rs10847864$ -e rs11610045$ /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > pos_controls.txt
cat LRRK2_HRC_coding_V2.txt pos_controls.txt > LRRK2_HRC_coding_V3.txt

#add in the first column ("ID") in CHR:POS format 
cut -f 1,2 LRRK2_HRC_coding_V3.txt | sed 's/\t/:/g' | sed 's/Chr:Start/ID/g' > first_column.txt
paste first_column.txt LRRK2_HRC_coding_V3.txt > LRRK2_HRC_coding_V4.txt

#copy the IDs back to working directory as LRRK2_coding_VOI.txt (VOI=variants of interest)
cut -f1 LRRK2_HRC_coding_V4.txt  | tail -n +2 > LRRK2_HRC_coding_V4_IDs.txt
scp LRRK2_HRC_coding_V4_IDs.txt /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt

# 48 variants present…
```

### 3.2 Perform GWAS of chromosome 12 for each cohort
note to self I will want to make this a batch script when I re-run

#### Normal GWAS for IPDGC cohorts
```
cd /data/LNG/Julie/Julie_LRRK2_Condi

mkdir NORMAL_GWAS_CHR12
module load plink/2.0-dev-20191128

cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt  | while read line
do 
	plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--chr 12 \
	--keep /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--pheno-name PHENO_PLINK --covar-variance-standardize \
	--pheno /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar-name AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.$line
done

## EXCEPTIONS => VANCE + MF no age...
plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.VANCE.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES//LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.VANCE

plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.MF

```

#### Conditional GWAS for IPDGC cohorts (no rs76904798 + no G2019S)
```
mkdir CONDI_GWAS_CHR12

cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
	plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--chr 12 \
	--keep /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.$line.txt \
	--pheno-name PHENOTYPE --covar-variance-standardize \
	--pheno /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.$line.txt \
	--covar /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.$line.txt \
	--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12.$line
done

## EXCEPTIONS => VANCE + MF no age...
plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.VANCE.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.VANCE.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.VANCE.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12.VANCE

plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.MF.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.MF.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_COVARIATES/LRRK2_condi_covariates.MF.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12.MF

```

#### Special conditional GWAS for IPDGC cohorts (no N2081D + no G2019S)
```
mkdir SPECIAL_GWAS_CHR12

cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
	plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--chr 12 \
	--keep /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.$line.txt \
	--pheno-name PHENOTYPE --covar-variance-standardize \
	--pheno /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.$line.txt \
	--covar /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.$line.txt \
	--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/SPECIAL_GWAS_CHR12.$line
done

## EXCEPTIONS => VANCE + MF no age...
plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.VANCE.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.VANCE.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.VANCE.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/SPECIAL_GWAS_CHR12.VANCE

plink2 --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--chr 12 \
--keep /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.MF.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.MF.txt \
--covar /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_COVARIATES/LRRK2_condi_covariates_SPECIAL.MF.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/SPECIAL_GWAS_CHR12.MF

```

### 3.3 Munging data (note this has been updated)

```
Files => 
NORMAL_GWAS_CHR12.*.PHENO_PLINK.glm.logistic.hybrid
SPECIAL_GWAS_CHR12.*.PHENOTYPE.glm.logistic.hybrid
CONDI_GWAS_CHR12.*.PHENOTYPE.glm.logistic.hybrid

#get rid of the first two columns (CHROM and POS) since they are repeats of ID
#keep ID, REF, ALT, A1, A1_FREQ, OR, LOG(OR)_SE, P
#columns of interest... => 3,4,5,6,13,19,20,22

### loop over to make the .hybrid files into readable .txt files 
cat /data/LNG/Julie/Julie_LRRK2_Condi/cohort_file.txt | while read line
do 
  cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/
  cut -f 3,4,5,6,13,19,20,22 NORMAL_GWAS_CHR12.$line.PHENO_PLINK.glm.logistic.hybrid > NORMAL_GWAS_CHR12.$line.txt
  cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12/
  cut -f 3,4,5,6,13,19,20,22 SPECIAL_GWAS_CHR12.$line.PHENOTYPE.glm.logistic.hybrid > SPECIAL_GWAS_CHR12.$line.txt
  cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/
  cut -f 3,4,5,6,13,19,20,22 CONDI_GWAS_CHR12.$line.PHENOTYPE.glm.logistic.hybrid > CONDI_GWAS_CHR12.$line.txt
done

#organize the files 
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
### 3.4 Pull LRRK2 coding variants from IPDGC CHR12 GWAS

Purpose: After incorporating UKB data, we will make forest plots of LRRK2 coding variants and positive controls.

```
##make new GWAS results files with only the LRRK2 coding VOI

cd /data/LNG/Julie/Julie_LRRK2_Condi/
cat cohort_file.txt | while read line
do 
  cd /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12
  #keep the header to add to VOI files
  head -1 NORMAL_GWAS_CHR12.$line.txt > header.txt
  #pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt NORMAL_GWAS_CHR12.$line.txt > temp.txt
  #combine these lines with the header
  cat header.txt temp.txt > NORMAL_GWAS_VOI.$line.txt
  rm header.txt temp.txt
  
  cd /data/LNG/Julie/Julie_LRRK2_Condi/SPECIAL_GWAS_CHR12
  #keep the header to add to VOI files
  head -1 SPECIAL_GWAS_CHR12.$line.txt > header.txt
  #pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt SPECIAL_GWAS_CHR12.$line.txt > temp.txt
  #combine these lines with the header
  cat header.txt temp.txt > SPECIAL_GWAS_VOI.$line.txt
  rm header.txt temp.txt
  
  cd /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12
  #keep the header to add to VOI files
  head -1 CONDI_GWAS_CHR12.$line.txt > header.txt
  #pull out the variants of interest
  grep -Ff /data/LNG/Julie/Julie_LRRK2_Condi/LRRK2_coding_VOI.txt CONDI_GWAS_CHR12.$line.txt > temp.txt
  #combine these lines with the header
  cat header.txt temp.txt > CONDI_GWAS_VOI.$line.txt
  rm header.txt temp.txt
done

##reorganize the files --> save the VOI files in LRRK2_coding_VOI within each GWAS directory

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
--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWASchr12.UKBB.EU.filtered_NEW \
--memory 235000

```
### 4.3 Subset phenotype data for covariate files

Coviariate files to be made: N=6
1) PD vs control 
2) PD vs control (no rs76904798 + no G2019S)
3) PD vs control (no N2081D + no G2019S)
4) PD proxy vs control 
5) PD proxy vs control (no rs76904798 + no G2019S)
6) PD proxy vs control (no N2081D + no G2019S)

Working directory: data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

#### Determine LRRK2 carrier status 

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS

module load plink/2.3-alpha
plink2 --bgen /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb_imp_chr12_v3.bgen --snps rs76904798,rs34637584,rs33995883 --make-bed --sample /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb33601_imp_chr1_v3_s487395.sample --out LRRK2_area_snps

module load plink
plink --bfile LRRK2_area_snps --recodeA --out LRRK2_area_snps2

LRRK2 carrier status stored here: LRRK2_area_snps2.raw

##Checking imputation quality:
#rs76904798_T => 5' risk variant
#rs34637584_A => G2019S 
#rs33995883_G => N2081D

grep 'rs76904798\|rs34637584\|rs33995883' /data/CARD/UKBIOBANK/IMPUTED_DATA/ukb_mfi_chr12_v3.txt

NAME		RS		BP		REF	ALT	MAF		ALT	R2		
12:40614434_C_T	rs76904798	40614434	C	T	0.144409	T	0.996952
rs34637584	rs34637584	40734202	G	A	0.000321086	A	1
rs33995883	rs33995883	40740686	A	G	0.0163887	G	1

```

#### Subset phenotype info based on LRRK2 status

```
module load R
cd data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
R
require("data.table")
require("dplyr")

### Read in all of the files:

#full covariate files
cov <- fread("/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt",header=T)
#PD cases 
PD <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD.txt",header=T)
#PD proxy cases
Proxy <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_parent_no_PD.txt",header=F)
#raw data filtered and unrelated -- use to filter the final groups
keep <- fread("/data/CARD/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.fam",header=F)
#Not controls: don't use PD cases or proxy cases as controls
not_control <- fread("/data/CARD/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_case_or_PD_parent.txt",header=T)
#Use LRRK2 status to subset covariate files for separate GWAS
LRRK2_status <- fread("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/LRRK2_area_snps2.raw",header=T)


#pull the sample IDs for each of the groups
PDshort <-  data.frame(PD$eid)
Proxyshort <- data.frame(Proxy$V1)
not_controlshort <- data.frame(not_control$eid)
keepshort <- data.frame(keep$V1)

#pull the relevant phenotype info
covshort <- data.frame(cov$FID, cov$AGE_OF_RECRUIT, cov$GENETIC_SEX)

#rename the sample IDs with "FID" 
names(PDshort) <- c("FID")
names(Proxyshort) <- c("FID")
names(not_controlshort) <- c("FID")
names(keepshort) <- c("FID")

#rename the cov file header
names(covshort) <- c("FID", "AGE", "SEX")

### Merge the phenotype info with the FIDs for each group to add age and sex info
PD2 = merge(PDshort,covshort,by.x='FID',by.y='FID')
Proxy2 = merge(Proxyshort,covshort,by.x='FID',by.y='FID')
not_control2 = merge(not_controlshort,covshort,by.x='FID',by.y='FID')

### Filter the groups using the filtered and unrelated data (keepshort)
PD3 = merge(PD2,keepshort,by.x='FID',by.y='FID')
Proxy3 = merge(Proxy2,keepshort,by.x='FID',by.y='FID')
not_control3 = merge(not_control2,keepshort,by.x='FID',by.y='FID')

### Get all controls....
keep3 = merge(keepshort,covshort,by.x='FID',by.y='FID')
#use anti_join to find unmatched records, those that are in keep3 but not not_control3 --> controls 
keep4 = anti_join(keep3,not_control3, by = c('FID'))
#we only want controls that are 60+ yrs old
keep5 = subset(keep4, AGE >= 60)

### Add column on PD status
PD3$STATUS <- "PD"
Proxy3$STATUS <- "PROXY"
keep5$STATUS <- "CONTROL"

### Subset controls randonly... 
#we want 10X as many controls as there are cases
# case-control = 1530 x 10 = 15300
nrow(PD3) #1530
#remember that keep5 is the controls >= 60 yrs
PD3_controls <- sample_n(keep5, 15300)
#use the rest of the controls for the Proxy GWAS
Proxy3_controls <- anti_join(keep5,PD3_controls, by = c('FID'))

### Cat dataframes...
PD_FINAL <- rbind(PD3, PD3_controls)
PROXY_FINAL <- rbind(Proxy3, Proxy3_controls)
#add the IID column
PD_FINAL$IID <- PD_FINAL$FID
PROXY_FINAL$IID <- PROXY_FINAL$FID 
#rearrange so that IID is the second column and age is last
PD_FINAL <- PD_FINAL[c(1,5,3,4,2)]
PROXY_FINAL <- PROXY_FINAL[c(1,5,3,4,2)]

### Subset based on LRRK2 carrier status
# rs76904798_T => 5' risk variant
# rs34637584_A => G2019S 
# rs33995883_G => N2081D

#get rid of unneeded columns: only need FID, rs76904798_T, rs34637584_A, rs33995883_G
LRRK2_status$IID <- NULL
LRRK2_status$PAT <- NULL
LRRK2_status$MAT <- NULL
LRRK2_status$SEX <- NULL
LRRK2_status$PHENOTYPE <- NULL
   
#add LRRK2 status to the final PD and Proxy datasets
PD_FINAL_LRRK2 <- merge(PD_FINAL,LRRK2_status,by.x='FID',by.y='FID')
PROXY_FINAL_LRRK2 <- merge(PROXY_FINAL,LRRK2_status,by.x='FID',by.y='FID')

#subset the final PD and Proxy datasets based on LRRK2  status
PD_FINAL_norisk_GS <- subset(PD_FINAL_LRRK2, rs76904798_T == 0 & rs34637584_A == 0)
PD_FINAL_noN2081D_GS <- subset(PD_FINAL_LRRK2, rs33995883_G == 0 & rs34637584_A == 0)
PROXY_FINAL_norisk_GS <- subset(PROXY_FINAL_LRRK2, rs76904798_T == 0 & rs34637584_A == 0)
PROXY_FINAL_noN2081D_GS <- subset(PROXY_FINAL_LRRK2, rs33995883_G == 0 & rs34637584_A == 0)

#save dataframes
write.table(PD_FINAL, file="UKB_PD_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_norisk_GS, file="UKB_PD_cases_control_over60_noriskGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_FINAL_noN2081D_GS, file="UKB_PD_cases_control_over60_noNDGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL, file="UKB_Proxy_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_norisk_GS, file="UKB_Proxy_cases_control_over60_noriskGS.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PROXY_FINAL_noN2081D_GS, file="UKB_Proxy_cases_control_over60_noNDGS.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

#organize some files
mkdir LRRK2_status_prep
mv LRRK2_area_snps* LRRK2_status_prep

```

### 4.4 Calculate PCs

See here for more information on flashpca: https://github.com/gabraham/flashpca

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
module load flashpca
module load plink

#test run flashpca on the PD normal group
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
#grep the subset filenames
ls UKB_P*> PC_files.txt

## make sure PC_files.txt contains these phenptype files: 
# UKB_PD_cases_control_over60_noNDGS.txt
# UKB_PD_cases_control_over60_noriskGS.txt
# UKB_PD_cases_control_over60.txt
# UKB_Proxy_cases_control_over60_noNDGS.txt
# UKB_Proxy_cases_control_over60_noriskGS.txt
# UKB_Proxy_cases_control_over60.txt

#use the phenotype files to subset each GWAS group for flashpca
#note this may take a few hours…
cat PC_files.txt  | while read line
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
# merge files in R
module load R
R
require(data.table)
require(dplyr)

#load the full covariates file
cov <- fread("/data/CARD/UKBIOBANK/ICD10_UKBB/Covariates/covariates_phenome_final.txt",header=T)
#remove the PC columns from cov so that we can merge with the new PCs
cov2 <- cov %>% select(1:8)

#pull the subset phenptype files from your working directory
file.names <- dir("/data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/", pattern="^UKB_P") 
for(file in file.names) {  
  pc <- fread(paste("pcs_",file,sep=""),header=T)
  pc$IID <- NULL
  pheno <- fread(file,header=T)
  #this info will become redundant when merging, keep (PD/CONTROL/PROXY) STATUS and LRRK2 carrier status
  pheno$IID <- pheno$SEX <- pheno$AGE <- NULL
  MM = merge(cov2,pheno,by='FID')
  MM2 = merge(MM,pc,by='FID')
  write.table(MM2, file=paste("COV_",file,sep=""), quote=FALSE,row.names=F,sep="\t")
}
q()
n


#replace PD/CONTROL/PROXY STATUS with numbers
cat /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/PC_files.txt | while read line
do 
  sed -i 's/PD/2/g' COV_$line
  sed -i 's/PROXY/2/g' COV_$line
  sed -i 's/CONTROL/1/g' COV_$line
done


#organize/remove some files
mkdir flashpca_files
mv eigenvalues_UKB_P* flashpca_files
mv eigenvectors_UKB_P* flashpca_files
mv pcs_UKB_P* flashpca_files
mv pve_UKB_P* flashpca_files
rm FILENAME_*
rm pruned_data*

```

### 4.5 Start GWAS on CHR 12 (note to self make the loop GWAS a bash script)

#### Perform GWAS in a loop

```
module load plink/2.0-dev-20191128

cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/UKB_GWAS/
mkdir GWAS_output
module load plink/2.0-dev-20191128

6 GWAS:
1) COV_UKB_PD_cases_control_over60_noNDGS.txt
2) COV_UKB_PD_cases_control_over60_noriskGS.txt
3) COV_UKB_PD_cases_control_over60.txt
4) COV_UKB_Proxy_cases_control_over60_noNDGS.txt
5) COV_UKB_Proxy_cases_control_over60_noriskGS.txt
6) COV_UKB_Proxy_cases_control_over60.txt

#test GWAS
# COV_UKB_PD_cases_control_over60_noNDGS.txt
# 1 binary phenotype loaded (1466 cases, 14775 controls).
plink2 --pfile chr12.UKBB.EU.filtered_NEW \
--pheno-name STATUS --pheno COV_UKB_PD_cases_control_over60_noNDGS.txt \
--covar COV_UKB_PD_cases_control_over60_noNDGS.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/COV_UKB_PD_cases_control_over60_noNDGS_chr12 --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize


###loop for all GWAS chr 12 only
#${line%%.*} allows you to remove the file extension
cat PC_files.txt  | while read line
do 
	plink2 --pfile chr12.UKBB.EU.filtered_NEW \
	--pheno-name STATUS --pheno COV_$line \
	--covar COV_$line --memory 235000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--out /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/COV_${line%%.*}_chr12 --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

```

#### Check the cases and controls included in each GWAS (note I want to ask about this)

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output/

### First check the case and control numbers in the covariate files 
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
[1] "COV_UKB_PD_cases_control_over60_noNDGS.txt has 1466 cases and 14839 controls"
[1] "COV_UKB_PD_cases_control_over60_noriskGS.txt has 1063 cases and 11232 controls"
[1] "COV_UKB_PD_cases_control_over60.txt has 1530 cases and 15300 controls"
[1] "COV_UKB_Proxy_cases_control_over60_noNDGS.txt has 12942 cases and 136250 controls"
[1] "COV_UKB_Proxy_cases_control_over60_noriskGS.txt has 9679 cases and 102931 controls"
[1] "COV_UKB_Proxy_cases_control_over60.txt has 13430 cases and 140908 controls"


### Compare this to the results from the GWAS log files
ls *.log > log_files.txt

cat log_files.txt  | while read line
do 
	echo $line
	grep "1 binary phenotype loaded" $line
done

Results: 
COV_UKB_PD_cases_control_over60_chr12.log
1 binary phenotype loaded (1529 cases, 15271 controls).
COV_UKB_PD_cases_control_over60_noNDGS_chr12.log
1 binary phenotype loaded (1466 cases, 14839 controls).
COV_UKB_PD_cases_control_over60_noriskGS_chr12.log
1 binary phenotype loaded (1063 cases, 11232 controls).
COV_UKB_Proxy_cases_control_over60_chr12.log
1 binary phenotype loaded (13404 cases, 140663 controls).
COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.log
1 binary phenotype loaded (12942 cases, 136250 controls).
COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.log
1 binary phenotype loaded (9679 cases, 102931 controls).

```

### 4.6 Make data ready for meta-analysis

These are the files to process:

COV_UKB_PD_cases_control_over60_chr12.STATUS.glm.logistic.hybrid
COV_UKB_PD_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid
COV_UKB_PD_cases_control_over60_noriskGS_chr12.STATUS.glm.logistic.hybrid
COV_UKB_Proxy_cases_control_over60_chr12.STATUS.glm.logistic.hybrid
COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid
COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.STATUS.glm.logistic.hybrid

#### Reformat plink2 GWAS output
```
#check the number of variants in these files…
wc -l COV_UKB_PD_cases_control_over60_noNDGS_chr12.STATUS.glm.logistic.hybrid 
# 1,357,686 variants--in comparison, 439,402 variants in /data/CARD/UKBIOBANK/FILTER_IMPUTED_DATA/chr12.UKBB.EU.filtered.pvar


cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/GWAS_output

ls COV_UKB_P*.hybrid | cat > GWAS_files.txt
module load python/3.6
cat GWAS_files.txt  | while read line
do 
	#filter the results by A1_FREQ (minor allele frequency) >=0.0001 --> to input.txt \
	awk '{ if($13 >= 0.0001) { print }}' $line > input.txt
	
	#reformat the plink2 results for meta-analysis using python
	if [[ $line == *"PD"* ]]; then
		python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
		--outfile toMeta.${line%%.*}.txt --B-or-C B
	elif [[ $line == *"Proxy"* ]]; then
		python /data/CARD/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile input.txt \
		--outfile toProxy.${line%%.*}.txt --B-or-C B; fi
done

#the PD files are ready for meta-analysis, the proxy files need more reformatting…

```

#### Adjust Proxy cases to the same scale as PD cases

```
# Make the toProxy files into .csv files
module load R
R
require("data.table")
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


# Convert the "normalized" proxy .csv files back to .txt files
module load R
R
require("data.table")
data1 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_chr12.csv",header=T)
data2 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.csv",header=T)
data3 <- fread("toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.csv",header=T)
write.table(data1, file="toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt",quote=F,row.names=F,sep="\t")
write.table(data2, file="toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt",quote=F,row.names=F,sep="\t")
write.table(data3, file="toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt",quote=F,row.names=F,sep="\t")
q()
n

```

#### Reformat UKB further and extract LRRK2 coding variants (note this is where I left off)

These are the files to reformat:

toMeta.COV_UKB_PD_cases_control_over60_chr12.txt
toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt
toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt
toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt
toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt
toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt

```
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS
mkdir METAL 
cd GWAS_output

#Copy over the reformatted toMeta files into the METAL directory
scp toMeta.*.txt /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/METAL

###we need the variant names of UKB to match those of IPDGC:
#EX of UKB case: 
head -2 /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/METAL/toMeta.COV_UKB_PD_cases_control_over60_chr12.txt 
markerID effectAllele alternateAllele beta se P effectAlleleFreq N rsID effectAlleleFreq_cases effectAlleleFreq_controls OR firthUsed error
12:87074:T:C T C 0.6356186982531077 0.793412 0.423062 0.00037558199999999996 16800 rs564020348 0.0006897289999999999 0.00034412900000000004 1.88819 N .
	
#Get rid of the :N:N off the end of the markerID for the UKB files…
cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/METAL
ls | grep UKB > list.txt
cat list.txt  | while read line
do
   	sed -i 's/:A//g' $line
	sed -i 's/:T//g' $line
	sed -i 's/:C//g' $line
	sed -i 's/:G//g' $line
done

# This is the format we need for making forest plots for all of the UKB toMeta files
# ID	REF	ALT	A1	A1_FREQ	OR	LOG(OR)_SE	P

cd /data/LNG/Julie/Julie_LRRK2_Condi/UKB_GWAS/METAL

```

#### Pull the LRRK2 coding variants from the UKB files for forest plots 

```
#note that the "special" conditional GWAS don't have the ND or GS variants

#create a .txt file with these variants and call it variants.txt

echo "12:40657700
12:40671989
12:40702911
12:40707778
12:40707861
12:40713899
12:40740686
12:40734202
12:40713901
12:40758652
12:40614434" > variants.txt



# UKB cases
head -1 toMeta.COV_UKB_PD_cases_control_over60_chr12.txt | cut -f 1-7 > header_PD.txt
grep -f variants.txt toMeta.COV_UKB_PD_cases_control_over60_chr12.txt | cut -f 1-7 > temp
cat header_PD.txt temp > NORMAL_GWAS_VOI.UKBPD.txt
grep -f variants.txt toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt | cut -f 1-7 > temp
cat header_PD.txt temp > CONDI_GWAS_VOI_SPECIAL.UKBPD.txt
grep -f variants.txt toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt | cut -f 1-7 > temp
cat header_PD.txt temp > CONDI_GWAS_VOI.UKBPD.txt


# UKB proxies
#we want to pull the b_adjusted, se_adjusted and p_derived for the proxy cases rather than b, se and p for the PD cases
head -1 toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt | cut -f 1,2,3,7,15,16,17 > header_proxy.txt
grep -f variants.txt toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt | cut -f 1,2,3,7,15,16,17 > temp
cat header_proxy.txt temp > NORMAL_GWAS_VOI.UKBproxy.txt
grep -f variants.txt toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt | cut -f 1,2,3,7,15,16,17 > temp
cat header_proxy.txt temp > CONDI_GWAS_VOI_SPECIAL.UKBproxy.txt
grep -f variants.txt toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt | cut -f 1,2,3,7,15,16,17 > temp
cat header_proxy.txt temp > CONDI_GWAS_VOI.UKBproxy.txt

```


#### ALSO SAVE THE FULL CHR12 GWAS RESULTS

```
###making the UKB CHR12 files…


# UKB cases
cut -f 1-7 toMeta.COV_UKB_PD_cases_control_over60_chr12.txt > NORMAL_GWAS_CHR12.UKBPD.txt
cut -f 1-7 toMeta.COV_UKB_PD_cases_control_over60_noNDGS_chr12.txt > CONDI_GWAS_CHR12_SPECIAL.UKBPD.txt
cut -f 1-7 toMeta.COV_UKB_PD_cases_control_over60_noriskGS_chr12.txt > CONDI_GWAS_CHR12.UKBPD.txt


# UKB proxies
#we want to pull the b_adjusted, se_adjusted and p_derived for the proxy cases rather than b, se and p for the PD cases
cut -f 1,2,3,7,15,16,17 toMeta.COV_UKB_Proxy_cases_control_over60_chr12.txt > NORMAL_GWAS_CHR12.UKBproxy.txt
cut -f 1,2,3,7,15,16,17 toMeta.COV_UKB_Proxy_cases_control_over60_noNDGS_chr12.txt > CONDI_GWAS_CHR12_SPECIAL.UKBproxy.txt
cut -f 1,2,3,7,15,16,17 toMeta.COV_UKB_Proxy_cases_control_over60_noriskGS_chr12.txt > CONDI_GWAS_CHR12.UKBproxy.txt


# copy to folders…
#note that I removed the dot before * so that SPECIAL files are sent to the CONDI folder
#copy the CHR12 files 
scp NORMAL_GWAS_CHR12* /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/
scp CONDI_GWAS_CHR12* /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/

#copy the VOI-filtered files
scp NORMAL_GWAS_VOI* /data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/VOI_files
scp CONDI_GWAS_VOI* /data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/VOI_files

```

#### UK Biobank files to combine with IPDGC for meta-analysis

```
...note that I need to change this section so that I put the SPECIAL UKB into the IPDGC SPECIAL GWAS folder

#all of CHR12 - proxy GWAS
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12_SPECIAL.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12.UKBproxy.txt

#all of CHR12 - PD GWAS
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/NORMAL_GWAS_CHR12.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12_SPECIAL.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/CONDI_GWAS_CHR12.UKBPD.txt

#variant of interest - proxy GWAS
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/VOI_files/NORMAL_GWAS_VOI.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/VOI_files/CONDI_GWAS_VOI_SPECIAL.UKBproxy.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/VOI_files/CONDI_GWAS_VOI.UKBproxy.txt

#variant of interest - PD GWAS
/data/LNG/Julie/Julie_LRRK2_Condi/NORMAL_GWAS_CHR12/VOI_files/NORMAL_GWAS_VOI.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/VOI_files/CONDI_GWAS_VOI_SPECIAL.UKBPD.txt
/data/LNG/Julie/Julie_LRRK2_Condi/CONDI_GWAS_CHR12/VOI_files/CONDI_GWAS_VOI.UKBPD.txt

```

## 5. Combining all data together

This section goes through: 
- Performing cohort level analyses on the created data from 1 and 4 and then meta-analyzing

``` 
# getting back to this:

cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/NORMAL/
HEADER of almost all:
ID	REF	ALT	A1	A1_FREQ	OR	LOG(OR)_SE	P


HEADER of UKB case:
markerID	effectAllele	alternateAllele	beta	se	P	effectAlleleFreq

HEADER of UKB proxy:
markerID	effectAllele	alternateAllele	effectAlleleFreq	b_adjusted	se_adjusted	p_derived

# remove column 2 from IPDGC and create beta
cat /data/LNG/CORNELIS_TEMP/LRRK2_conditional/cohort_file.txt | while read line
do 
	Rscript --vanilla quick_reformat.R $line
done

# reformat UKB
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' NORMAL_GWAS.UKBPD.txt > NORMAL_GWAS.UKBPDv2.txt
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' NORMAL_GWAS.UKBproxy.txt > NORMAL_GWAS.UKBproxyv2.txt

# create files per sample

head -1 NORMAL_GWAS.DUTCHv2.txt > header.txt
cat variants.txt  | while read line
do 
   	grep $line NORMAL_GWAS.*v2.txt > $line.txt
	cat header.txt $line.txt > header_$line.txt
	sed -e 's/NORMAL_GWAS.//g' header_$line.txt | sed -e 's/.txt:'$line'//g' > header_"$line"v2.txt
done

```
Now doing the same for => CONDITIONAL

```
cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/CONDI/
HEADER of almost all:
ID	REF	ALT	A1	A1_FREQ	OR	LOG(OR)_SE	P


HEADER of UKB case:
markerID	effectAllele	alternateAllele	beta	se	P	effectAlleleFreq

HEADER of UKB proxy:
markerID	effectAllele	alternateAllele	effectAlleleFreq	b_adjusted	se_adjusted	p_derived

# remove column 2 from IPDGC and create beta
cat /data/LNG/CORNELIS_TEMP/LRRK2_conditional/cohort_file.txt | while read line
do 
	Rscript --vanilla quick_reformat.R $line
done

# reformat UKB
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $7, $4, $5, $6}' CONDI_GWAS.UKBPD.txt > CONDI_GWAS.UKBPDv2.txt
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $2, $4, $5, $6, $7}' CONDI_GWAS.UKBproxy.txt > CONDI_GWAS.UKBproxyv2.txt

# create files per sample

head -1 CONDI_GWAS.DUTCHv2.txt > header.txt
cat variants.txt  | while read line
do 
   	grep $line CONDI_GWAS.*v2.txt > $line.txt
	cat header.txt $line.txt > header_$line.txt
	sed -e 's/CONDI_GWAS.//g' header_$line.txt | sed -e 's/.txt:'$line'//g' > header_"$line"v2.txt
done

```


###### making forest plots...


```
POS:change
Rscript --vanilla ../forest_plot_LRRK2.R 12:40657700 p.Asn551Lys
Rscript --vanilla ../forest_plot_LRRK2.R 12:40671989 p.Ile723Val
Rscript --vanilla ../forest_plot_LRRK2.R 12:40702911 p.Arg1398His
Rscript --vanilla ../forest_plot_LRRK2.R 12:40707778 p.Arg1514Gln
Rscript --vanilla ../forest_plot_LRRK2.R 12:40707861 p.Pro1542Ser
Rscript --vanilla ../forest_plot_LRRK2.R 12:40713899 p.Met1646Thr
Rscript --vanilla ../forest_plot_LRRK2.R 12:40740686 p.Asn2081Asp
Rscript --vanilla ../forest_plot_LRRK2.R 12:40734202 p.Gly2019Ser
Rscript --vanilla ../forest_plot_LRRK2.R 12:40713901 p.Ser1647Thr
Rscript --vanilla ../forest_plot_LRRK2.R 12:40758652 p.Met2397Thr
Rscript --vanilla ../forest_plot_LRRK2.R 12:40614434 rs76904798

Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40657700 p.Asn551Lys
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40671989 p.Ile723Val
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40702911 p.Arg1398His
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40707778 p.Arg1514Gln
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40707861 p.Pro1542Ser
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40713899 p.Met1646Thr
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40740686 p.Asn2081Asp
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40734202 p.Gly2019Ser
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40713901 p.Ser1647Thr
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40758652 p.Met2397Thr
Rscript --vanilla ../forest_plot_LRRK2_condi.R 12:40614434 rs76904798


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# start like this
# Rscript --vanilla forest_plot_LRRK2.R $FILENAME $FILENAME2
# Rscript --vanilla forest_plot_LRRK2.R 12:40713899 Met1646Thr
FILENAME = args[1]
FILENAME2 = args[2]
print(args[1])
print(args[2])
print(FILENAME)
print(FILENAME2)
library(metafor)
data <- read.table(paste("header_",FILENAME,"v2.txt",sep=""), header = T)
##data <- read.table("header_12:40713899v2.txt", header = T)
labs <- data$ID
yi   <- data$beta
sei  <- data$LOG.OR._SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file = paste(FILENAME,"_final.pdf",sep=""), width = 8, height = 6)
Pvalue <- formatC(resFe$pval, digits=4, format="f")
## pdf(file = "12:40713899_final.pdf", width = 8, height = 6)
forest(resFe, xlim=c(-2,2), main=paste(FILENAME2," P=",Pvalue ,sep=""),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.5,0.75, 1, 2, 3)))
dev.off()
#png(file = paste(FILENAME,"_final.png",sep=""), width = 8, height = 6)
#forest(resFe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.9931472), main=paste(FILENAME2," option",sep=""),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.5,0.75, 1, 2, 3)))
#dev.off()

```

### done... for now ....


#### CHECK SNCA variant to confirm things....
4:90641340 => rs356220

```
module load plink/2.0-dev-20191128

cd /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/


cat /data/LNG/CORNELIS_TEMP/LRRK2_conditional/cohort_file.txt  | while read line
do 
	plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--snp 4:90641340 \
	--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--pheno-name PHENO_PLINK --covar-variance-standardize \
	--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.$line.txt \
	--covar-name AGE,SEX_COV,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/NORMAL_GWAS.$line
done

## EXCEPTIONS => VANCE + MF no age...

plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snp 4:90641340 \
--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.VANCE.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.VANCE.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/NORMAL_GWAS.VANCE

plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snp 4:90641340 \
--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--pheno-name PHENO_PLINK --covar-variance-standardize \
--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/COVARIATES/LRRK2_condi_covariates_NORMAL.MF.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/NORMAL_GWAS.MF


### condi
cat /data/LNG/CORNELIS_TEMP/LRRK2_conditional/cohort_file.txt | while read line
do 
	plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
	--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
	--snp 4:90641340 \
	--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.$line.txt \
	--pheno-name PHENOTYPE --covar-variance-standardize \
	--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.$line.txt \
	--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.$line.txt \
	--covar-name AGE,SEX,PC1,PC2,PC3,PC4,PC5 \
	--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/CONDI_GWAS.$line
done

## EXCEPTIONS => VANCE + MF no age...

plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snp 4:90641340 \
--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.VANCE.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.VANCE.txt \
--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.VANCE.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/CONDI_GWAS.VANCE

plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --memory 99000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--snp 4:90641340 \
--keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.MF.txt \
--pheno-name PHENOTYPE --covar-variance-standardize \
--pheno /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.MF.txt \
--covar /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_covariates.MF.txt \
--covar-name SEX,PC1,PC2,PC3,PC4,PC5 \
--out /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/CONDI_GWAS.MF

##
cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/
grep 90641340 CONDI_GWAS.*hybrid > CONDI_signal.txt
grep 90641340 NORMAL_GWAS.*hybrid > NORMAL_signal.txt
head -1 CONDI_GWAS.VANCE.PHENOTYPE.glm.logistic.hybrid > header.txt

```

```
Make forest plots
cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/GWAS_NEW/SNCA/

CONDI_forest.txt
NORMAL_forest.txt

module load R
R
library(metafor)
data <- read.table("CONDI_forest.txt", header = T)
## data <- read.table("NORMAL_forest.txt", header = T) 
labs <- data$DATA
yi   <- data$BETA
sei  <- data$LOG.OR._SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file = "SNCA_Final_condi.pdf", width = 8, height = 6)
# pdf(file = "SNCA_Final_normal.pdf", width = 8, height = 6)
Pvalue <- formatC(resFe$pval, digits=4)
forest(resFe, xlim=c(-2,2), main=paste("SNCA variant P=",Pvalue ,sep=""),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.5,0.75, 1, 2, 3)))
dev.off()
#png(file = paste(FILENAME,"_final.png",sep=""), width = 8, height = 6)
#forest(resFe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.9931472), main=paste(FILENAME2," option",sep=""),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9, at=log(c(0.5,0.75, 1, 2, 3)))
#dev.off()



```





## 6. Check LD co-inheritance of LRRK2 coding variants
This section goes through:
- First check which coding variants can be imputed here...?
- Checking if LRRK2 variants are typically co-inherited or not?
- LRRK2 G2019S with all other coding variants
- LRRK2 rs76904798 with all other coding variants
- Preparing files for Tables for manuscript


### 6.1 - First check which coding variants can be imputed here...?

```
mkdir HRC_LRRK2
head -1 /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > header.txt
awk '$1 == "12"' /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt | grep LRRK2 | grep exonic | sed 's/ /_/g' | grep nonsynonymous > temp
cat header.txt temp > LRRK2_HRC_coding.txt
# manually add in 12:40614434 / rs76904798
grep rs76904798 /data/CARD/GENERAL/HRC_ouput_annovar_ALL.txt > risk_variant.txt
cat LRRK2_HRC_coding.txt risk_variant.txt > LRRK2_HRC_coding_V2.txt
cut -f 1,2 LRRK2_HRC_coding_V2.txt | sed 's/\t/:/g' | sed 's/Chr:Start/ID/g' > first_row.txt
paste first_row.txt LRRK2_HRC_coding_V2.txt > LRRK2_HRC_coding_V3.txt
# 45 variants present...

```

### 6.2 - Assessing frequency of these variants in full data-set...

```
# Create a couple subset files to use....
# G2019S carriers only
plink --bfile HARDCALLS_PD_september_2018_no_cousins --snps 12:40734202 --recodeA --out LRRK2_G2019S_only

module load R
R
data <- read.table("LRRK2_G2019S_only.raw",header=T)
newdata <- subset(data, X12.40734202_A == 0) 
dim(newdata) # 33757     7
# adding some additional sampleinfo
cov <- read.table("IPDGC_all_samples_covariates.txt",header=T)
# drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL
MM = merge(newdata,cov,by='FID')
dim(MM) # 33757    43
write.table(MM, file="LRRK2_G2019S_only_with_COV.txt", quote=FALSE,row.names=F,sep="\t")

# rs76904798 carriers only
plink --bfile HARDCALLS_PD_september_2018_no_cousins --snps 12:40614434 --recodeA --out LRRK2_rs76904798_only
module load R
R
data <- read.table("LRRK2_rs76904798_only.raw",header=T)
newdata <- subset(data, X12.40614434_T == 0) 
dim(newdata) # 33335     7
# adding some additional sampleinfo
cov <- read.table("IPDGC_all_samples_covariates.txt",header=T)
# drop some columns because otherwise merge conflict
cov$IID <- NULL
cov$fatid <- NULL
cov$matid <- NULL
MM = merge(newdata,cov,by='FID')
dim(MM) # 33335    43
write.table(MM, file="LRRK2_rs76904798_only_with_COV.txt", quote=FALSE,row.names=F,sep="\t")

# Ecluding all G2019S AND rs76904798 carriers
LRRK2_condi_sample_selection.txt

# copy back to workingfolder
scp LRRK2_G2019S_only_with_COV.txt /data/LNG/CORNELIS_TEMP/LRRK2_conditional/HRC_LRRK2/
scp LRRK2_rs76904798_only_with_COV.txt /data/LNG/CORNELIS_TEMP/LRRK2_conditional/HRC_LRRK2/

```

```
module load plink
# all data 
# 21478 are cases and 24388 are controls.
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --assoc --out freq
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --logistic --out freq
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --model --out freq

# excluding LRRK2 G2019S carriers and NA
# 16072 are cases and 17685 are controls.
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --assoc --out freq_no_G2019S --keep LRRK2_G2019S_only_with_COV.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --logistic --out freq_no_G2019S --keep LRRK2_G2019S_only_with_COV.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --model --out freq_no_G2019S --keep LRRK2_G2019S_only_with_COV.txt

# excluding '5 risk variant carriers and NA
# 15318 are cases and 18017 are controls.
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --assoc --out freq_no_rs76904798 --keep LRRK2_rs76904798_only_with_COV.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --logistic --out freq_no_rs76904798 --keep LRRK2_rs76904798_only_with_COV.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --model --out freq_no_rs76904798 --keep LRRK2_rs76904798_only_with_COV.txt

# excluding both...
# 11441 are cases and 13091 are controls.
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --assoc --out freq_no_G2019S_rs76904798 --keep LRRK2_condi_sample_selection.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --logistic --out freq_no_G2019S_rs76904798 --keep LRRK2_condi_sample_selection.txt
plink --bfile /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/HARDCALLS_PD_september_2018_no_cousins \
--extract LRRK2_HRC_coding_V3.txt --model --out freq_no_G2019S_rs76904798 --keep LRRK2_condi_sample_selection.txt
```

So many files!! now merge these in R

```
module load R
R
# normal
data1 <- read.table("freq.assoc",header=T)
data2 <- read.table("freq.assoc.logistic",header=T)
data3 <- read.table("freq.model",header=T)
MM <- data1[,c(2,4,5,6,7)]
MM2 <- data2[,c(2,6,7,9)]
MM3 <- subset(data3, TEST=="GENO")
MM4 <- MM3[,c(2,6,7)]
MM5 <- merge(MM,MM4,by="SNP")
normal <- merge(MM5,MM2, by="SNP")
colnames(normal) <- c("SNP","A1","F_A_normal","F_U_normal","A2","AFF_normal","UNAFF_normal","NMISS_normal","OR_normal","P_normal")
write.table(normal, file="LRRK2_coding_variants_freq_normal.txt", quote=FALSE,row.names=F,sep="\t")

# no G2019S
data1 <- read.table("freq_no_G2019S.assoc",header=T)
data2 <- read.table("freq_no_G2019S.assoc.logistic",header=T)
data3 <- read.table("freq_no_G2019S.model",header=T)
MM <- data1[,c(2,4,5,6,7)]
MM2 <- data2[,c(2,6,7,9)]
MM3 <- subset(data3, TEST=="GENO")
MM4 <- MM3[,c(2,6,7)]
MM5 <- merge(MM,MM4,by="SNP")
noGS <- merge(MM5,MM2, by="SNP")
colnames(noGS) <- c("SNP","A1","F_A_noGS","F_U_noGS","A2","AFF_noGS","UNAFF_noGS","NMISS_noGS","OR_noGS","P_noGS")
write.table(noGS, file="LRRK2_coding_variants_freq_noGS.txt", quote=FALSE,row.names=F,sep="\t")

# no rs76904798
data1 <- read.table("freq_no_rs76904798.assoc",header=T)
data2 <- read.table("freq_no_rs76904798.assoc.logistic",header=T)
data3 <- read.table("freq_no_rs76904798.model",header=T)
MM <- data1[,c(2,4,5,6,7)]
MM2 <- data2[,c(2,6,7,9)]
MM3 <- subset(data3, TEST=="GENO")
MM4 <- MM3[,c(2,6,7)]
MM5 <- merge(MM,MM4,by="SNP")
noRisk <- merge(MM5,MM2, by="SNP")
colnames(noRisk) <- c("SNP","A1","F_A_noRisk","F_U_noRisk","A2","AFF_noRisk","UNAFF_noRisk","NMISS_noRisk","OR_noRisk","P_noRisk")
write.table(noRisk, file="LRRK2_coding_variants_freq_noRisk.txt", quote=FALSE,row.names=F,sep="\t")

# no G2019S and rs76904798
data1 <- read.table("freq_no_G2019S_rs76904798.assoc",header=T)
data2 <- read.table("freq_no_G2019S_rs76904798.assoc.logistic",header=T)
data3 <- read.table("freq_no_G2019S_rs76904798.model",header=T)
MM <- data1[,c(2,4,5,6,7)]
MM2 <- data2[,c(2,6,7,9)]
MM3 <- subset(data3, TEST=="GENO")
MM4 <- MM3[,c(2,6,7)]
MM5 <- merge(MM,MM4,by="SNP")
noGSRisk <- merge(MM5,MM2, by="SNP")
colnames(noGSRisk) <- c("SNP","A1","F_A_noGSRisk","F_U_noGSRisk","A2","AFF_noGSRisk","UNAFF_noGSRisk","NMISS_noGSRisk","OR_noGSRisk","P_noGSRisk")
write.table(noGSRisk, file="LRRK2_coding_variants_freq_noGSRisk.txt", quote=FALSE,row.names=F,sep="\t")

```

Then inspect closer and assess whether there are differences....

Two variants most affected by removal of risk variants are:

LRRK2:NM_198578:exon42:c.A6241G:p.N2081D

LRRK2:NM_198578:exon32:c.G4541A:p.R1514Q

Table with frequencies below...		

| SNP | 12:40740686/N2081D | 12:40707778/R1514Q |  NOTE |
| ------------- | ------------- | ------------- | ------------- |
| Freq_PD_normal | 0.02361 | 0.009172 | |
| Freq_Control_normal | 0.0189 | 0.008037 | |
| Freq_PD_noGS | 0.02535 | 0.009115 | |
| Freq_Control_noGS | 0.02075 | 0.007181 |
| Freq_PD_noRisk | 0.003101 | 0.0008813 | Big Drop... |
| Freq_Control_noRisk | 0.002525 | 0.001499 | Big Drop... |
| Freq_PD_noGSRisk | 0.00354 | 0.0007429 | Big Drop... |
| Freq_Control_noGSRisk	| 0.002559 | 0.001222 | Big Drop... |


### 6.3 - Preparing files for Tables for manuscript

Lets first do this for the IPDGC data....

```
# IPDGC data (FULL)
cd /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/
# note using hardcall imputed data R2>0.8
module load plink
# create variant of interest file
echo 12:40657700 > LRRK2_variants_Nicole.txt
echo 12:40671989 >> LRRK2_variants_Nicole.txt
echo 12:40702911 >> LRRK2_variants_Nicole.txt
echo 12:40707778 >> LRRK2_variants_Nicole.txt
echo 12:40707861 >> LRRK2_variants_Nicole.txt
echo 12:40713899 >> LRRK2_variants_Nicole.txt
echo 12:40734202 >> LRRK2_variants_Nicole.txt
echo 12:40740686 >> LRRK2_variants_Nicole.txt
echo 12:40713901 >> LRRK2_variants_Nicole.txt
echo 12:40758652 >> LRRK2_variants_Nicole.txt
echo 12:40614434 >> LRRK2_variants_Nicole.txt
# check model output
plink --bfile HARDCALLS_PD_september_2018_no_cousins --model --extract LRRK2_variants_Nicole.txt --out LRRK2_variants_Nicole_model
# extract GENO from model
grep GENO LRRK2_variants_Nicole_model.model
# check MAF per disease status
plink --bfile HARDCALLS_PD_september_2018_no_cousins --assoc --extract LRRK2_variants_Nicole.txt --out LRRK2_variants_Nicole_assoc
# results including all data (no variant exclusion)
Variant		Cases		Controls		F_A		F_U
12:40614434	558/5602/15318	486/5885/18017		0.1564		0.1406
12:40657700	93/2588/18797	139/3110/21139		0.06458		0.06946
12:40671989	130/3062/18286	140/3230/21018		0.07733		0.07196
12:40702911	95/2642/18741	142/3152/21094		0.06593		0.07044
12:40707778	3/388/21087	1/390/23997		0.009172	0.008037
12:40707861	15/1167/20296	31/1439/22918		0.02787		0.03077
12:40713899	12/866/20600	5/820/23563		0.02072		0.01702
12:40713901	1775/8885/10818	2250/10067/12071	0.2895		0.2987
12:40734202	1/236/16072	0/20/17685		0.007297	0.0005648
12:40740686	11/992/20475	13/896/23479		0.02361		0.0189
12:40758652	2561/9757/9160	2817/10786/10785	0.3464		0.3366

# IPDGC data (no GS and no '5 carriers)
plink --bfile HARDCALLS_PD_september_2018_no_cousins --model --extract LRRK2_variants_Nicole.txt --out LRRK2_variants_Nicole_model_no5_noGS --keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/HRC_LRRK2/LRRK2_condi_sample_selection.txt
# extract GENO from model
grep GENO LRRK2_variants_Nicole_model_no5_noGS.model
# check MAF per disease status
plink --bfile HARDCALLS_PD_september_2018_no_cousins --assoc --extract LRRK2_variants_Nicole.txt --out LRRK2_variants_Nicole_assoc_no5_noGS --keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/HRC_LRRK2/LRRK2_condi_sample_selection.txt
# results excluding GS and '5 risk variant
Variant		Cases		Controls		F_A		F_U
12:40614434	0/0/11441	0/0/13091		0		0
12:40657700	67/1615/9759	95/1945/11051		0.07644		0.08154
12:40671989	107/1924/9410	91/2007/10993		0.09344		0.08361
12:40702911	67/1640/9734	98/1969/11024		0.07753		0.08269
12:40707778	0/17/11424	0/32/13059		0.0007429	0.001222
12:40707861	10/678/10753	13/828/12250		0.0305		0.03262
12:40713899	8/550/10883	4/518/12569		0.02474		0.02009
12:40713901	1228/5039/5174	1584/5751/5756		0.3276		0.3407
12:40734202	0/0/11441	0/0/13091		0		0
12:40740686	0/81/11360	0/67/13024		0.00354		0.002559
12:40758652	1602/5432/4407	1694/5907/5490		0.3774		0.355
```

```
Another way of doing this

echo 12:40657700 > variants.txt
echo 12:40671989 >> variants.txt
echo 12:40702911 >> variants.txt
echo 12:40707778 >> variants.txt
echo 12:40707861 >> variants.txt
echo 12:40713899 >> variants.txt       
echo 12:40740686 >> variants.txt      
echo 12:40734202 >> variants.txt
echo 12:40713901 >> variants.txt
echo 12:40758652 >> variants.txt
echo 12:40614434 >> variants.txt

cd /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/

plink --bfile HARDCALLS_PD_september_2018_no_cousins --extract variants.txt --freq

plink --bfile HARDCALLS_PD_september_2018_no_cousins --extract variants.txt \
--assoc --keep /data/LNG/CORNELIS_TEMP/LRRK2_conditional/LRRK2_condi_sample_selection.txt


```


OK next lets first do this for the UK Biobank exome data....

```
# IPDGC data (FULL)
cd /data/LNG/CORNELIS_TEMP/LRRK2_conditional/UKB/
module load plink
# create variant of interest file
echo 12:40263898:C:G > LRRK2_variants_Nicole.txt
echo 12:40278187:A:G >> LRRK2_variants_Nicole.txt
echo 12:40309109:G:A >> LRRK2_variants_Nicole.txt
echo 12:40313976:G:A >> LRRK2_variants_Nicole.txt
echo 12:40314059:C:T >> LRRK2_variants_Nicole.txt
echo 12:40320097:T:C >> LRRK2_variants_Nicole.txt
echo 12:40340400:G:A >> LRRK2_variants_Nicole.txt
echo 12:40346884:A:G >> LRRK2_variants_Nicole.txt
# using LRRK2_ALL_variants_euro_subset_mac10.* from below...
plink --bfile LRRK2_ALL_variants_euro_subset_mac10 --extract LRRK2_variants_Nicole.txt --freq --out LRRK2_variants_Nicole_freq
# results europeans not excluding anyone
SNP		MAF	NCHROBS
12:40263898:C:G	0.06485	82494
12:40278187:A:G	0.06611	82498
12:40309109:G:A	0.06775	82498
12:40313976:G:A	0.008885	82498
12:40314059:C:T	0.03159	82498
12:40320097:T:C	0.01855	82498
12:40340400:G:A	0.0003152	82498
12:40346884:A:G	0.01484	82498
# now check same data excluding '5 carriers and GS carriers
module load R
R
data <- read.table("LRRK2_ALL_variants_euro_subset_mac10_RECODE.raw",header=T)
dim(data) # 41249    43
# '5 risk => rs76904798_T
# GS => X12.40340400.G.A_A
newdata <- subset(data, X12.40340400.G.A_A == 0) 
dim(newdata) # 41223     43
newdata2 <- subset(newdata, rs76904798_T == 0) 
dim(newdata2) # 29970    43
NEW <- data.frame(newdata2$FID,newdata2$IID,newdata2$SEX)
colnames(NEW) <- c("FID","IID","SEX")
write.table(NEW, file="LRRK2_exome_euro_no5andGS_carriers.txt", quote=FALSE,row.names=F,sep="\t")
# check freq in plink
plink --bfile LRRK2_ALL_variants_euro_subset_mac10 --extract LRRK2_variants_Nicole.txt --freq --out LRRK2_variants_Nicole_freq_no5GS --keep LRRK2_exome_euro_no5andGS_carriers.txt
# results europeans excluding '5 and G2019S carriers
SNP		MAF	NCHROBS
12:40263898:C:G	0.07431	59936
12:40278187:A:G	0.07608	59940
12:40309109:G:A	0.07754	59940
12:40313976:G:A	0.002386	59940
12:40314059:C:T	0.0373	59940
12:40320097:T:C	0.02137	59940
12:40340400:G:A	0	59940
12:40346884:A:G	0.001835	59940
```
Done...


## 7 - Checking R1441C and R1441G carrier status in AMP-PD....

```
cd /PATH/TO/PD/AMP-PD/VCFs/

chr12:40310434 C  T 	rs33939927 R-C
chr12:40310434 C  G	rs33939927 R-G
https://www.snpedia.com/index.php/Rs33939927

module load vcftools
module load bcftools
# extract only R1441 regions from VCF file
vcftools --gzvcf chr12.vcf.gz --out LRRK2_R1441_status --chr chr12 --from-bp 40310333 --to-bp 40310535 --recode --recode-INFO-all
# split multiallelics option -m add in - for splitting and add in + for joining
bcftools norm -m- LRRK2_R1441_status.recode.vcf > LRRK2_R1441_status_split_allele.recode.vcf
```

## Done....


![myImage](https://media.giphy.com/media/XRB1uf2F9bGOA/giphy.gif)




