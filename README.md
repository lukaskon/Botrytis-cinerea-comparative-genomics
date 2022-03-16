# BotrytisDNA_pipeline


# FastQC
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=60G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name runFastQC     # you can give your job a name for easier identification (same as -J)
#SBATCH -o FastQC_trimmed_slurm2

########## Command Lines to Run ##########

module load fastqc                   ### load necessary modules, e.g.

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_DNA
mkdir FastQC_MultiQC_reports


fastqc *fastq.gz -o FastQC_MultiQC_reports/



scontrol show job $SLURM_JOB_ID 
```

# MultiQC
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=60G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name runMultiQC     # you can give your job a name for easier identification (same as -J)


########## Command Lines to Run ##########


module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load MultiQC                   ### load necessary modules, e.g.

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisRNASeq/RawReadsRNA/20211109_mRNASeq_PE150/MultiQC

multiqc .


scontrol show job $SLURM_JOB_ID
```

# Fastp (initial trim)
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=36:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=8G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name trim_fastp      # you can give your job a name for easier identification (same as -J)
#SBATCH -o Fastp_DNA_slurm

########## Command Lines to Run ##########

conda install -c bioconda fastp

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/RawReadsDNA/20211029_DNASeq_PE150

mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_DNA/Fastp_reports

for infile in *_R1_001.fastq.gz

do

base=$(basename ${infile} _R1_001.fastq.gz)
fastp -i ${infile} -I ${base}_R2_001.fastq.gz \
-o /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_DNA/${base}_R1_trim_UP.fastq.gz -O /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_D
NA/${base}_R2_trim_UP.fastq.gz \
--adapter_fasta adapters.fa -f 15 -F 15 -l 36 -q 20 \
-h /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_DNA/Fastp_reports/${base}_fastp.html -j /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trim
med_DNA/Fastp_reports/${base}_fastp.json

done

scontrol show job $SHOW_JOB_ID     ### write job information to output file



########## Comments #################


#For adapters, refer to https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa
```




# FastQC and MultiQC again


see above scripts


# Obtain reference genome


http://ftp.ensemblgenomes.org/pub/fungi/release-52/fasta/botrytis_cinerea/dna/


Bcinerea_RefGenome_Ensembl$ ls
Botrytis_cinerea.ASM83294v1.dna.toplevel.fa          Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.ann
Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.fai      Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.bwt.2bit.64
Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.0123  Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.bwt.8bit.32
Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.amb   Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz.pac


# Index reference genome with Bwa-Mem2
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=05:00:00 # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1 # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=30G # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Bcinerea_index # you can give your job a name for easier identification (same as -J)
#SBATCH -o Index_Bcinerea_slurm

########## Command Lines to Run ##########

module load icc/2018.1.163-GCC-6.4.0-2.28  impi/2018.1.163
module load bwa-mem2/2.0                   ### load necessary modules, e.g.

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/DNAscripts/Bcinerea_RefGenome_Ensembl

/opt/software/bwa-mem2/2.0pre2/bwa-mem2 index Botrytis_cinerea.ASM83294v1.dna.toplevel.fa.gz


scontrol show job $SLURM_JOB_ID
```


# Add read group information here?



# Align samples to reference genome with Bwa-Mem2
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=30G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Aln_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o Align_Bcinerea_slurm

########## Command Lines to Run ##########

module load icc/2018.1.163-GCC-6.4.0-2.28  impi/2018.1.163
module load bwa-mem2/2.0                   ### load necessary modules, e.g.

mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/BwaMem2_aligned_Bcinerea

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Fastp_Trimmed_DNA




for infile in *_R1_trim_UP.fastq.gz

do

base=$(basename ${infile} _R1_trim_UP.fastq.gz)

/opt/software/bwa-mem2/2.0pre2/bwa-mem2 mem -t 32 /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/DNAscripts/Bcinerea_RefGenome_Ensembl/Botrytis_cinerea.ASM83294v1.dna
.toplevel.fa.gz ${base}_R1_trim_UP.fastq.gz ${base}_R2_trim_UP.fastq.gz > /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/BwaMem2_aligned_Bcinerea/${base}_aln_bwamem2.
sam


done

scontrol show job $SLURM_JOB_ID 
```


# Convert Sam to Bam
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=8           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name SamToBam_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o SamToBam_slurm

########## Command Lines to Run ##########

module load GCC/9.3.0
module load SAMtools/1.11


mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Bam_BwaMem2_aligned_Bcinerea
cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/BwaMem2_aligned_Bcinerea  ### change to the directory where your data is located



for infile in *_L002_aln_bwamem2.sam

do

base=$(basename ${infile} _L002_aln_bwamem2.sam)

samtools view --threads 8 -Sb ${base}_L002_aln_bwamem2.sam > /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Bam_BwaMem2_aligned_Bcinerea/${base}_aln_bwamem2.bam

done

scontrol show job $SLURM_JOB_ID 
```

# Sort and index Bam files
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########



#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Ce_samToBam      # you can give your job a name for easier identification (same as -J)
#SBATCH -o SortIndexBam_slurm



########## Command Lines to Run ##########



module load GCC/9.3.0
module load SAMtools/1.11

mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/SortInx_Bam_Bcinerea
cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Bam_BwaMem2_aligned_Bcinerea




for infile in *.bam

do

base=$(basename ${infile} .bam)
samtools sort --threads 32 ${base}.bam -o ${base}_sort.bam
samtools index ${base}_sort.bam

mv ${base}_sort.bam /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/SortInx_Bam_Bcinerea

done


scontrol show job $SLURM_JOB_ID
```

# Mark and remove duplicates
```
#!/bin/bash --login
########### Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=80G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Picard_dups      # you can give your job a name for easier identification (same as -J)
#SBATCH -o Picard_removedups_slurm

########### Command Lines to Run ##########


module load picard/2.18.1-Java-1.8.0_152


cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/5_SortInx_Bam_Bcinerea

for infile in *.bam

do

base=$(basename ${infile} .bam)
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${base}.bam O=/mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea/${base}_rmvdups.bam M=/
mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea/metrics/${base}_rmvdups_metrics.txt ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=
50000 REMOVE_DUPLICATES=true

done

scontrol show job $SLURM_JOB_ID 
```

# Use GATK IndelRealigner here?

# Calculating coverage across all samples

more runSamtools_cov.sb
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=8           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name SamCov_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o SamCoverage_slurm

########## Command Lines to Run ##########

module load GCC/9.3.0
module load SAMtools/1.11


mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea/coverage
cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea


for infile in *.bam

do

base=$(basename ${infile} .bam)

samtools coverage ${base}.bam -o /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea/coverage/${base}.txt

done

scontrol show job $SLURM_JOB_ID     ### write job information to output file
```



command line: 

coverage$ awk 'FNR>1{a+=$6;b++}END{print "Average: " a/b}' *txt
Average: 95.3505

Using samtools depth -a

depth$ awk 'FNR>1{a+=$3;b++}END{print "Average: " a/b}' AF15_S92_aln_bwamem2_sort_rmvdups_depth.txt
Average: 21.6842



For all 96: ranges from 11-49x, with an average of 24x coverage depth.




# Call variants

```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=20G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Bcftools_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o Variant_calls_correct_slurm

########## Command Lines to Run ##########

module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load bcftools/1.9.64

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/6_Duplicates_removed_Bcinerea

bcftools mpileup -a AD,DP,SP -Ou -f /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/0_DNAscripts/Bcinerea_RefGenome_Ensembl/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa *.bam|bcftools call -f GQ,GP -mO z -o /mn
t/research/Hausbeck_group/Lukasko/BotrytisDNASeq/unfilteredVCF/Bcin_unfiltered.vcf.gz


scontrol show job $SLURM_JOB_ID     ### write job information to output file 
```


# Take a random sample of variants and run stats to decide on filtering parameters

### To view the number of variants in a file:
```
module load bcftools
bcftools view -H I16_S58_aln_bwamem2_sort_rmvdups_calls.vcf.gz |wc -l
```

## Pull a random 200,000 variants, then compress and index

```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=40G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name subsetVar      # you can give your job a name for easier identification (same as -J)
#SBATCH -o UnfilteredVCF_subset_index_slurm

########## Command Lines to Run ##########


cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/7_unfilteredVCF
mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/8_unfilteredVCF_subset


module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load bcftools/1.9.64
conda install -c bioconda vcflib



bcftools view Bcin_unfiltered.vcf.gz |vcfrandomsample -r 0.0046 > /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/8_VCF_subset_Bcinerea_subset/Bcin_unfiltered_subset.vcf

# sample rate was chosen as 0.0024 so it would yeild ~200,000 random VCFs (200,000/42,801,380=~0.0046)


echo Random 200,000 variants pulled.

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/8_unfilteredVCF_subset

#compress (like bgzip, not gzip)
bcftools view Bcin_unfiltered_subset.vcf -Oz -o Bcin_unfiltered_subset.vcf
echo Compressed.

#index compressed files to make it easier to access for running stats later
bcftools index Bcin_unfiltered_subset.vcf.gz
echo Indexed.


# https://github.com/vcflib/vcflib/blob/master/doc/vcfrandomsample.md
# https://speciationgenomics.github.io/filtering_vcfs/

scontrol show job $SLURM_JOB_ID     ### write job information to output file
```


## Stats
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=18:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=32           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=60G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name Vcftools_Bcin      # you can give your job a name for easier identification (same as -J)
#SBATCH -o VCFstats_filtered_slurm

########## Command Lines to Run ##########

module load GNU/7.3.0-2.30  OpenMPI/3.1.1
module load VCFtools/0.1.15-Perl-5.28.0


cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/7_unfilteredVCF
mkdir /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats

vcftools --gzvcf Bcin_filtered.vcf.gz --freq2 --out /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz --max-alleles 2

vcftools --gzvcf Bcin_filtered.vcf.gz --depth --out /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz

vcftools --gzvcf Bcin_filtered.vcf.gz --site-mean-depth --out  /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz

vcftools --gzvcf Bcin_filtered.vcf.gz --site-quality --out /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz

vcftools --gzvcf Bcin_filtered.vcf.gz --missing-indv --out /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz

vcftools --gzvcf Bcin_filtered.vcf.gz --missing-site --out /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/9_FilteredVCF_Stats/Bcin_filteredMAF.vcf.gz


# MAF in output file name indicates that this was filtered with MAF of 0.1


scontrol show job $SLURM_JOB_ID     ### write job information to output file
```

## Move stats output onto local machine and use R to view. (See Variant Based Statistics.R in current repository)


