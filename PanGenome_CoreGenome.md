# Creating fastas of 286 Botrytis cinerea isolates

## Use vcf to create fasta
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########


#SBATCH --time=12:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=50G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=makeFASTAs      # you can give your job a name for easier identification (same as -J)
#SBATCH -o CreateFASTAs1_slurm

########## Command Lines to Run ##########

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Extractions/vcf

for infile in *.vcf
do

base=$(basename ${infile} .vcf)

module load GCCcore/6.4.0
module load tabix/0.2.6

bgzip -c ${base}.vcf > ../../FASTAs/${base}.vcf.gz
cd ../../FASTAs
tabix ${base}.vcf.gz

module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load bcftools/1.9.64

bcftools consensus -f /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/0_DNAscripts/ReferenceGenome/Botrytis_cinerea.ASM83294v1.dna.toplevel.
fa ${base}.vcf.gz > ${base}.fa

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/Extractions/vcf

done

scontrol show job $SHOW_JOB_ID

```

## Annotate proteins with Augustus (in batches because it takes approx 1hr/fasta)
```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########


#SBATCH --time=100:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=80G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=augustus_annotation_qz      # you can give your job a name for easier identification (same as -J)
#SBATCH -o Augustus_QZ_slurm

########## Command Lines to Run ##########

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/FASTAs

module purge
conda activate funannotate

#module load GCC/11.3.0  OpenMPI/4.1.4
#module load AUGUSTUS/3.5.0
#AUGUSTUS_CONFIG_PATH=/opt/software/AUGUSTUS/3.5.0-foss-2022a/config


for infile in [Q-Z]*.fa

do

base=$(basename ${infile} .fa)

augustus ${base}.fa --species=botrytis_cinerea > augustus_annotations/${base}_Augustus.gff

done

scontrol show job $SLURM_JOB_ID

```
for file in *; do echo ${file} >> ProteinCounts_augustus.txt; grep -rn "protein" ${file}|wc -l >> ProteinCounts_augustus.txt; done

## Identify orthologous group with OrthoMCL

#### First, request database and MySQL access: https://docs.icer.msu.edu/MySQL_configuration/

orthomcl.config provided by HPCC:
```

```
https://docs.icer.msu.edu/Load_the_software/ 
https://docs.icer.msu.edu/orthomcl-pipeline/

```
module load icc/2016.3.210-GCC-5.4.0-2.26  impi/5.1.3.181
module load OrthoMCL/2.0.9-Perl-5.24.0

orthomclPairs orthomcl.config log_file cleanup=[yes|no|only|all] <startAfter=TAG>
```










