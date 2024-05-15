
## Running fastStructure using structure_threader
https://structure-threader.readthedocs.io/en/latest/

*~/fastStructure* - path to fastStructure binary file


### Conversion from VCF to PED
```
vcftools --vcf Carper.vcf --plink --out test_plink
```
### Conversion from PED to BED
```
plink --file test_plink --make-bed --out test_bed
```
### Running structure threader 

```
structure_threader run -K 10 -R 1 -i test_bed.bed -o test_output/ -t 4 -fs ~/fastStructure --ind popmap.txt
```

*selecting only biallelic loci:*
 *awk '/^#/ || $5 !~ /,/' test.vcf > test.2alleles.vcf*  
