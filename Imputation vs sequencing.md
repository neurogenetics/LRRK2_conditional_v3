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



```
