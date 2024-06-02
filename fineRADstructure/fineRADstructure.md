

### Download fineRADstructure-tools
```
git clone https://github.com/edgardomortiz/fineRADstructure-tools
```


### Convert file from ipyrad format 
```
python fineRADstructure/finerad_input.py --input FILENAME.alleles.loci
```

### Calculate the co-ancestry matrix: 
```
RADpainter paint FILENAME.alleles.loci.min2.finerad
```



### Assign individuals to populations: 
```
~/fineRADstructure/finestructure -x 100000 -y 100000 -z 1000 FILENAME.alleles.loci.min2_chunks.out FILENAME_chunks.mcmc.xml
```


### Tree building: 
```
~/fineRADstructure/finestructure -m T -x 10000 FILENAME.alleles.loci.min2_chunks.out FILENAME_chunks.mcmc.xml FILENAME_chunks.mcmcTree.xml
```

  

Then plot in R using FinestructureLibrary.R and fineRADstructurePlot.R (https://github.com/millanek/fineRADstructure)
