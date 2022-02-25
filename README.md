# BotrytisDNA_pipeline


# FastQC


# MultiQC


# Fastp (initial trim)


# FastQC and MultiQC again


# Obtain reference genome


# Index reference genome with Bwa-Mem2


# Align samples to reference genome with Bwa-Mem2



# Convert Sam to Bam, then sort and index

# ...other steps


# Calculating coverage across all samples

more runSamtools_cov.sb
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



on command line: 
coverage$ awk 'FNR>1{a+=$6;b++}END{print "Average: " a/b}' *txt
Average: 95.3505



