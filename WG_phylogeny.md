# Create a maximum-likelihood phylogenetic tree using WGS from 276 isolates + reference

#### Start with fastas that were created in "PanGenome_CoreGenome.md"

### Compile and index

cat *.fa > AllFastas_withRef.fasta
samtools faidx AllFastas_withRef.fasta

### Multiple alignment with MAFFT

```
module load GCC/10.2.0  OpenMPI/4.0.5
module load MAFFT/7.475-with-extensions

mafft --auto --thread 16 AllFastas_withRef.fasta > AllFastas_aligned_test.fasta
```
