# Plink: filtering VCF and preparing for data visualization
==============================================================================
### Refer to "Statistical Population Genomics" 
### Chapter 4: Data Management and Summary Statistics with PLINK by Christopher C. Chang


Note: Everything in this markdown file can be done command line quickly.

### Validate variants on command line (takes a ~5m)
```
BotrytisDNASeq$ gatk-4.2.5.0/gatk ValidateVariants -V 10_FilteredVCF/Plates123/PLINK/BcinereaP123.SNVonly.filteredPASS_renamed.vcf
```

Ensure that VCF contains mono or bi-allelic only

```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=03:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=30G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Bcftools_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o BCF_biallelic_slurm

########## Command Lines to Run ##########

module load bcftools/1.9.64

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/10_FilteredVCF/Plates123/PLINK/


bcftools view --max-alleles 2 --exclude-types indels BcinereaP123.SNVonly.filteredPASS_renamed.vcf > Bcin123_biallelic.vcf



scontrol show job $SLURM_JOB_ID     ### write job information to output file
```

load plink 1.9

```
module spider plink #version needed will depend on task
```

## Convert VCF to BED file (plink binary format) 


```C++
plink2 --vcf BcinereaP123.SNVonly.filteredPASS_renamed.vcf \
--make-bed \
--out Bcin123_SNV1
```


## Rename variants without a unique ID for easy downstream processing

```
plink --bfile Bcin123_SNV1 \
--set-missing-var-ids @:# \
--make-bed \
--out Bcin123_SNV_idfilled2
```

No pedigree information to add.


## Filtering for missingness

```
plink --bfile Bcin123_SNV_idfilled2 \
--geno 0.1 \ #throws out variants where >10% of the calls are "NA"s
--mind 0.1 \ #throws out every samples where >10% of the genotype calls are "NA"s
--make-bed \
--out Bcin123_SNV_MissFiltered3
```

## Minor allele freq reporting and filtering out variants with MAF < 5%

```
plink --bfile Bcin123_SNV_MissFiltered3 \
--freq \
--out Bcin123_allelefreq

plink --bfile Bcin123_SNV_MissFiltered3 \
--maf 0.05 \
--make-bed \
--out Bcin123_MAF05
  
```
711121 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
570854 variants and 278 people pass filters and QC.

## Check for deviation from Hardy-Weinberg equilibrium
```
plink2 --bfile Bcin123_MAF05 \
--hwe 1e-25 keep-fewhet \
--make-bed \
--out Bcin123_hwe4
```
--hwe keep-fewhet: 0 variants removed due to Hardy-Weinberg exact test
(founders only).
570854 variants remaining after main filters.


## Selecting a SNP subet in approc linkage equilibrium

```
plink2 --bfile Bcin123_hwe4 \
--indep-pairwise 200kb 1 0.5 \
--out Bcin123_LDpruned
```
--indep-pairwise (16 compute threads): 538735/570854 variants removed.

## Create SNP subset
```
plink --bfile Bcin123_hwe4 --extract Bcin123_LDpruned.prune.in --make-bed --out Bcin123_LDsubset
```
32119 variants and 278 people pass filters and QC.

## PCA
```
plink2 --bfile Bcin123_LDsubset --pca 5 --out BcinMAF05_PCA_results
```

Download locally to view in R


```
setwd("C:/Users/nikki/Michigan State University/PSM.Hausbecklab - Nikki Lukasko - Nikki Lukasko/Botrytis/Molecular/Population Structure Analysis/Plink PCA")

pca_table <- read.table("BcinMAF05_PCA_results.eigenvec", header = TRUE, comment.char = "")
plot(pca_table[, c("PC1", "PC2", "PC3", "PC4", "PC5")])
```

## Export to VCF for other processes (may want LD pruned or one step before)

```
plink2 --bfile Bcin123_LDsubset --ref-from-fa /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/0_DNAscripts/ReferenceGenome/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa --export vcf --out Bcin123_plinkLDpruned
```


