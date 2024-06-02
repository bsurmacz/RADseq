### Change name of the file (!important - new extension of files in ipyrad -> .loci instead of alleles.loci): 

```
mv FILENAME.loci FILENAME.alleles.loci 
```

### Download fineRADstructure-tools
```
git clone https://github.com/edgardomortiz/fineRADstructure-tools
```


### Convert file from ipyrad format 
```
python fineRADstructure-tools/finerad_input.py --input FILENAME.alleles.loci
```

### Calculate the co-ancestry matrix: 
```
RADpainter paint FILENAME.alleles.loci.min2.finerad
```



### Assign individuals to populations: 
```
finestructure -x 100000 -y 100000 -z 1000 FILENAME.alleles.loci.min2_chunks.out FILENAME_chunks.mcmc.xml
```


### Tree building: 
```
finestructure -m T -x 10000 FILENAME.alleles.loci.min2_chunks.out FILENAME_chunks.mcmc.xml FILENAME_chunks.mcmcTree.xml
```

  

Then plot in R using FinestructureLibrary.R and fineRADstructurePlot.R (https://github.com/millanek/fineRADstructure)
