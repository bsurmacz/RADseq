
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


### NOW FOR ALL FILES
```
conda activate vcftools
mkdir PEDs
for file in biallelic/*.vcf; do
    base=$(basename "$file" .vcf)

    # Run vcftools and send output to PEDs directory
    vcftools --vcf "$file" --plink --out "PEDs/${base}"
done
```

```
mkdir BEDs
conda activate plink
for file in PEDs/*.ped; do
    base=$(basename "$file" .ped)
    plink --file PEDs/"$base" --make-bed --out "BEDs/${base}"
done
```

##
## : generate popmaps!!
```
mkdir BEDs

conda activate bcftools
for file in biallelic/*.vcf; do
    base=$(basename "$file" .ped)
    bcftools query -l "$file" | awk -F'-' '{print $0, $1}' > "$base".popmap
done



```

```
for file in BEDs/*.bed; do
    base=$(basename "$file" .bed)
   structure_threader run -K 20 -R 1 -i BEDs/"$base".bed -o structure_output/"$base" -t 16 -fs /opt/miniconda3/bin/fastStructure --ind "$base".vcf.popmap
done
```

 
