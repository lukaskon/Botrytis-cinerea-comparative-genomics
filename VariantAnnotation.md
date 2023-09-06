

# snpEff variant annotation


Genome file is part of snpeff config file. Run command line.
```
java -Xmx4g -jar snpEff/snpEff.jar Botrytis_cinerea 10_FilteredVCF/Plates123/BcinereaP123.SNV_INDEL.filteredPASS_renamed.vcf > 10_FilteredVCF/Plates123/BcinereaP123.SNV_INDEL.filteredPASS_renamed_ann.vcf
```


### Run with group S and sensu scricto separately
subset with:
```
bcftools view --samples-file GroupS/Pop2_GroupS --min-ac=1 --no-update BcinereaP123.SNV_INDEL.filteredPASS_renamed.vcf > GroupSsubset.vcf
```
annotate with:
```
java -Xmx8g -jar snpEff/snpEff.jar Botrytis_cinerea -stats 11_Finals/SNPEFF/GroupS/ -csvStats 11_Finals/SNPEFF/GroupS/GroupSsubset.vcf > 11_Finals/SNPEFF/GroupS/GroupSsubset_ann.vcf
```
