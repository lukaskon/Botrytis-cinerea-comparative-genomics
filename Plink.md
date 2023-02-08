
# Plink: filtering VCF and preparing for data visualization
==============================================================================


### Validate variants on command line (takes a ~5m)
```
BotrytisDNASeq$ gatk-4.2.5.0/gatk ValidateVariants -V 10_FilteredVCF/Plates123/PLINK/BcinereaP123.SNVonly.filteredPASS_renamed.vcf
```
```
module spider plink #version needed will depend on task
```

## Convert VCF to BED file (plink binary format) 


```C++
plink2 --vcf BcinereaP123.SNVonly.filteredPASS_renamed.vcf \
--make-bed \
--out Bcin123_SNV1.bed
```


## Rename variants without a unique ID for easy downstream processing

```
plink --bfile Bcin123_SNV1.bed \
--set-missing-var-ids @:# \
--make-bed \
--out Bcin123_SNV_idfilled2.bed
```

No pedigree information to add.


## Filtering for missingness

```
plink --bfile Bcin123_SNV_idfilled2.bed
--geno 0.1 \ #throws out variants where >10% of the calls are "NA"s
--mind 0.1 \ #throws out every samples where >10% of the genotype calls are "NA"s
--make-bed \
--out Bcin123_SNV_MissFiltered3.bed
```

## Minor allele freq reporting and filtering out variants with MAF < 5%

```
plink -bfile Bcin123_SNV_MissFiltered3.bed
--freq \
--out Bcin123_allelefreq

plink -bfile Bcin123_SNV_MissFiltered3.bed
--maf 0.05 \
--make-bed \
--out Bcin123_MAF05.bed
  
```












