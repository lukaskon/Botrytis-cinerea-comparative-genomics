# ADMIXTURE analysis


Starting with .bed file and other outputs for reading from plink filtering (Bcin123_LDsubset*)

```
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=60G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name runADMIXTURE_SNP     # you can give your job a name for easier identification (same as -J)
#SBATCH -o AdmixtureSNPs_chooseK_slurm

########## Command Lines to Run ##########

module load admixture                   ### load necessary modules, e.g.

cd /mnt/research/Hausbeck_group/Lukasko/BotrytisDNASeq/PopStrucutre_admixture/Plates123


for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv -j16 Bcin123_LDsubset.bed $K | tee log${K}.out; done

scontrol show job $SLURM_JOB_ID

```

```
grep -h CV log*.out
```
Choose the number (K) that has the lowest cross-validation error value.





