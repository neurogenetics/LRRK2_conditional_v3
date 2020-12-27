## Checking how good is the imputation vs sequencing data

```
# Julie and Cornelis
# December 2020
# Two datasets
1) UK Biobank => genotypes vs exome seq
2) IPDGC data and LNG WGS => genotypes vs genome seq

working dir:
cd /data/CARD/projects/LRRK2_sandbox/

```

### UK Biobank

UK Biobank => genotypes vs exome seq


##### Processing imputed genotype data (hg19)
```
# using imputed genotyped from:
cd /data/CARD/UKBIOBANK/IMPUTED_DATA/

# subset LRRK2 area
module load plink/2.3-alpha

plink2 --bgen ukb_imp_chr12_v3.bgen --chr 12 --from-kb 40604 --to-kb 40769 \
--make-pgen --out LRRK2_area --sample ukb33601_imp_chr1_v3_s487395.sample 
# 5224 variants remaining after main filters.

# copy to working dir
scp LRRK2_area.* /data/CARD/projects/LRRK2_sandbox

# working in sandbox
cd /data/CARD/projects/LRRK2_sandbox/

plink2 --pfile LRRK2_area --make-bed --out LRRK2_area_imputed

plink2 --bfile LRRK2_area_imputed --extract RS_to_keep.txt --make-bed --out LRRK2_area_imputed_RS_only
# 10 variants remaining after main filters.

# for merge and keep OG name
plink2 --bfile LRRK2_area_imputed --extract RS_to_keep.txt --make-bed --out LRRK2_area_imputed_for_merge --keep LRRK2_area_exome_for_merge.fam

```

##### Processing exome sequencing data
```
# using exome sequencing data from:
cd /data/CARD/UKBIOBANK/EXOME_DATA_200K/PLINK_files/

module load plink/2.3-alpha
plink2 --bed ukb23155_c12_b0_v1.bed --fam ukb23155_c1_b0_v1_s200632.fam \
--bim UKBexomeOQFE_chr12.bim --chr 12 --from-kb 40210 --to-kb 40375 --out LRRK2_area_exome --make-bed
# 3692 variants remaining after main filters.

# copy to working dir
scp LRRK2_area_exome.* /data/CARD/projects/LRRK2_sandbox

# working in sandbox
cd /data/CARD/projects/LRRK2_sandbox/

grep -f extract_exome_variants.txt LRRK2_area_exome.bim | cut -f 2 > exome_to_RS.txt

# add manually rs numbers to exome_to_RS.txt + remove 12:40309109:G:T and 12:40346884:A:T

plink2 --bfile LRRK2_area_exome --extract exome_to_RS.txt --update-name exome_to_RS.txt --make-bed --out LRRK2_area_exome_RS_only
# 10 variants remaining after main filters.

# for merge and keep OG name
plink2 --bfile LRRK2_area_exome --extract exome_to_RS.txt --make-bed --out LRRK2_area_exome_for_merge

```

##### Combining data
```
module load plink #1.9
# Report all mismatching calls.
plink --bfile LRRK2_area_exome_RS_only --bmerge LRRK2_area_imputed_RS_only --merge-mode 6
# 10 markers
# 2000800 overlapping calls, 1997001 nonmissing in both filesets.
# 1994359 concordant, for a concordance rate of 0.998677.

# deeper dive...
plink --bfile LRRK2_area_exome_for_merge --bmerge LRRK2_area_imputed_for_merge --recodeA --out MERGED_DATA

```

##### Results table

numbers from excel but Im sure there is a quicker way of doing this...

total => the sum of 0, 1, 2 from MERGED_DATA.raw => rs numbers 

confirmed => concatenation of imputed (0 ,1 , 2) and exome (0 ,1 , 2) and count of 00, 11 and 22


| RS-ID  | AA change | total/confirmed REF 0 | total/confirmed HET 1 | total/confirmed ALT 2 |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| rs33995463  | L119P | 199057/198998 | 1011/822 | 1/1 | 
| rs7308720  | N551K | 173697/173653 | 25448/25180 | 924/918 | 
| rs10878307  | I723V | 174384/174384 | 24700/24681 | 985/981 | 
| rs7133914  | R1398H | 172695/172681 | 26410/25906 | 964/963 | 
| rs35507033  | R1514Q | 196515/196508 | 3534/3490 | 20/20 | 
| rs33958906  | P1542S | 187626/187527 | 10895/10832 | 158/158 | 
| rs35303786  | M1646T | 191948/191920 | 6068/5942 | 52/49 | 
| rs11564148  | S1647T | 98418/98249 | 83780/83150 | 17871/17798 | 
| rs33995883  | N2081D | 193490/193486 | 6517/6150 | 62/62 | 
| rs3761863  | M2397T | 88603/88496 | 89043/88858 | 22423/22415 | 



### IPDGC data

IPDGC data and LNG WGS => genotypes vs genome seq

```
# done previously a while ago...
check => OVERLAP_GENOMES_AND_IPDGC_GENOTYPE_DATA.xlsx

ID stored:
WGS => WGS_overlap_with_IPDGC.txt
IPDGC => IPDGC_overlap_with_WGS.txt

```

##### Processing imputed genotype data (hg19)
```
# using imputed genotyped from:
cd /data/LNG/CORNELIS_TEMP/PD_FINAL_PLINK_2018/
HARDCALLS_PD_september_2018_no_cousins.

module load plink/2.3-alpha

plink2 --bfile HARDCALLS_PD_september_2018_no_cousins --keep /data/CARD/projects/LRRK2_sandbox/IPDGC_overlap_with_WGS.txt \
--make-bed --out /data/CARD/projects/LRRK2_sandbox/IPDGC_overlap_WGS --mac 1 
# 1947 samples

# subset variants of interest
cd /data/CARD/projects/LRRK2_sandbox/
plink2 --bfile IPDGC_overlap_WGS --extract IPDGC_variants_to_RS.txt --update-ids update_fam_IPDGC.txt --keep WGS_overlap_IPDGC.fam \
--out IPDGC_for_merge --make-bed 
# 1899 samples remaining

plink2 --bfile IPDGC_for_merge --extract IPDGC_variants_to_RS.txt --update-name IPDGC_variants_to_RS.txt --keep WGS_overlap_IPDGC.fam \
--out IPDGC_for_merge_check --make-bed 

```

##### Processing genome sequencing data
```
# using genome sequencing data from:
cd /data/CARD/PD/GENOMES/august19/genotypes/

module load plink/2.3-alpha
plink2 --pfile pd.june2019.chr12.sqc --keep /data/CARD/projects/LRRK2_sandbox/WGS_overlap_with_IPDGC.txt \
--make-bed --out /data/CARD/projects/LRRK2_sandbox/WGS_overlap_IPDGC --mac 1 
# 1899 samples

# subset variants of interest
cd /data/CARD/projects/LRRK2_sandbox/
plink2 --bfile WGS_overlap_IPDGC --extract WGS_variants_to_RS.txt  \
--out WGS_for_merge --make-bed 
# 1899 samples remaining

plink2 --bfile WGS_overlap_IPDGC --extract WGS_variants_to_RS.txt --update-name WGS_variants_to_RS.txt \
--out WGS_for_merge_check --make-bed 

```
##### Combining data
```
module load plink #1.9
# Report all mismatching calls.
plink --bfile IPDGC_for_merge_check --bmerge WGS_for_merge_check --merge-mode 6
# 10 markers
# 18990 overlapping calls, 18972 nonmissing in both filesets.
# 18966 concordant, for a concordance rate of 0.999684.

# deeper dive...
plink --bfile IPDGC_for_merge --bmerge WGS_for_merge --recodeA --out MERGED_DATA_IPDGC

```

##### Results table

numbers from excel but Im sure there is a quicker way of doing this...

total => the sum of 0, 1, 2 from MERGED_DATA_IPDGC.raw => rs numbers 

confirmed => concatenation of imputed (0 ,1 , 2) and WGS (0 ,1 , 2) and count of 00, 11 and 22


| RS-ID  | AA change | total/confirmed REF 0 | total/confirmed HET 1 | total/confirmed ALT 2 |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| rs33995463  | L119P | 1873/1873 | 8/8 | 0/0 | 
| rs7308720  | N551K | 1664/1664 | 226/225 | 9/9 | 
| rs10878307  | I723V | 1599/1598 | 290/290 | 10/10| 
| rs7133914  | R1398H | 1656/1656 | 233/233 | 10/10 | 
| rs35507033  | R1514Q | 1861/1861 | 38/37 | 0/0 | 
| rs33958906  | P1542S | 1792/1792 | 106/106 | 1/1 | 
| rs35303786  | M1646T | 1831/1829 | 68/68 | 0/0 | 
| rs11564148  | S1647T | 906/906 | 830/829 | 163/163 | 
| rs33995883  | N2081D | 1816/1816 | 82/82 | 1/1 | 
| rs3761863  | M2397T | 840/840 | 837/837 | 222/222 | 

## Done

Main conclusion => Imputation is pretty accurate in comparison with WGS and exome seq

