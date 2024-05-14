library(vcfR)
library(ggplot2)
library(adegenet)
library(poppr)
library(dplyr)
library(hierfstat)
library(vegan)
library(leaflet)

#################################################
###  Definitions of functions:
################################################


# filter loci from genind object, returns filtered genind
filter_loci<-function(gen,METADATA, MAF_filter=0.0,
                     max_heterozygosity_per_locus_filter=1,
                     missing_per_locus_filter=1,
                     missing_per_individual_filter=1,
                     number_of_alleles_per_locus=2){
  
  #distribution of number of alleles per SNP
  PLOT_histogram_number_of_alleles<-  ggplot()+aes(as.vector(nAll(gen)))+geom_histogram()+ylab("number of loci")+xlab("Alleles per locus")

  #keeping loci with no more than 2 alleles
  loci_to_keep <- nAll(gen) <= number_of_alleles_per_locus
  gen <- gen[loc = loci_to_keep]
  
  # maf filter
  gen<-informloci(gen,cutoff=0, MAF = MAF_filter)
  
  # missing data per locus filter
  gen<-missingno(gen, type = "loci", cutoff = missing_per_locus_filter, quiet = FALSE, freq = FALSE)
  
  #missing data per individual filter
  gen<-missingno(gen, type = "geno", cutoff = missing_per_individual_filter, quiet = FALSE, freq = FALSE)
  
  #SAM
  samples_populations<-data.frame(Sample=rownames(gen$tab))
  samples_populations<-samples_populations %>% left_join(METADATA,by="Sample" )
  # conversion from genind to hierfstat
  Data.hierfstat<-genind2hierfstat(gen,pop=factor(samples_populations$Population))
  
  STATS<-basic.stats(Data.hierfstat)
  Het<-STATS$perloc
  
  
  # histogram of heterozygosity per locus
  PLOT_histogram_heterozygosity_per_locus<-  ggplot()+aes(Het$Ho)+geom_histogram()+xlab("Heterozygosity per locus")
  
  # maximum heterozygosity per locus filter 
  het_filtered <- Het[which(Het$Ho  < max_heterozygosity_per_locus_filter), ]
  het_loci <- as.vector(rownames(het_filtered))
  het_loci <- gsub('X', '', het_loci)
  gen.het <- gen[loc=het_loci]
  
  return(list(gen.het,PLOT_histogram_number_of_alleles,PLOT_histogram_heterozygosity_per_locus))
  
  }

population_statistics<-function(gen.het,METADATA){
  
  samples_populations<-data.frame(Sample=rownames(gen.het$tab))
  samples_populations<-samples_populations %>% left_join(METADATA,by="Sample" )

    # converting to hierfstat
  Data.hierfstat<-genind2hierfstat(gen.het,pop=factor(samples_populations$Population))
  STATS<-basic.stats(Data.hierfstat)
  print(STATS) 
  # mean observed heterozygosity 
  Ho_populations<-apply(STATS$Ho,2,function(x){mean(x,na.rm=T)})
  Hs_populations<-apply(STATS$Hs,2,function(x){mean(x,na.rm=T)})
  Fis_populations<-apply(STATS$Fis,2,function(x){mean(x,na.rm=T)})
  # plotting PCA of allele frequencies in individuals
  PCA<-indpca(Data.hierfstat)
  plot(PCA)
  
  #genetic distance matrix
  Fst<-genet.dist(Data.hierfstat) #,method="WC84")  #
  hist(Fst)
  
  # Multi Dimensional Scaling plot of distances among populations
  pcplotM<-metaMDS(Fst)
  plot(pcplotM,"sites")
  orditorp(pcplotM,"sites",cex=1.8)
  
  # Writting genetic distance matrix to file
  #genetic_distance_matrix<-as.matrix(Fst)+t(as.matrix(Fst))
  
  
  populations_info<-METADATA %>% group_by(Population) %>% summarise(Latitude=mean(Latitude),Longitude=mean(Longitude))
  MDS_POINTS<-data.frame(pcplotM$points)
  MDS_POINTS$Population<-rownames(MDS_POINTS)
  pop_results<-data.frame(Ho=Ho_populations,Hs=Hs_populations,Fis=Fis_populations,Population=names(Fis_populations))
  pop_results<-pop_results %>% left_join(MDS_POINTS, by="Population")
  Populations_results_merged<-populations_info %>% left_join(pop_results, by="Population")
  Populations_results_merged<-Populations_results_merged %>% filter(!is.na(Ho))
  
  return(Populations_results_merged)
  }


plot_map<-function(data,var){
  pal<-colorRampPalette(c("red","green"))
  pal<-colorNumeric(palette = pal(100), domain = data[var])
  
  map_MDS <- leaflet(data = data) %>%
    # add a tile layer
    addTiles() %>%
    # add circle markers based on latitude and longitude columns, colored by "Value"
    addCircleMarkers(
      lng = ~Longitude, lat = ~Latitude,
      fillColor = ~pal( unlist(data[var])),fillOpacity =1,
      label = ~as.character(paste("Population",`Population`,var,": ", unlist(data[var]) )),
      radius=10) %>% addLegend("bottomleft", pal = pal, values = ~unlist(data[var]))
  map_MDS
}


##########################################################
##########################################################
#### Script starts:

# input data:
vcf.test<-"test.vcf"
sample_data<-read.csv("test_sample_data.csv")


# reading and converting
Test_vcf<-read.vcfR(vcf.test)
Test_gen<-vcfR2genind(Test_vcf)

# filtering data 
Test_gen_filtered<-filter_loci(Test_gen,sample_data,MAF_filter = 0.0,
                               missing_per_locus_filter=0.6,
                               missing_per_individual_filter=0.7)
# generating basic stats
Test_STATS<-population_statistics(Test_gen_filtered[[1]],sample_data)

# N alleles per locus:
Test_gen_filtered[[2]]
# Heterozygosity per locus:
Test_gen_filtered[[3]]


# Showing coordinates on MDS axes on map:
plot_map(Test_STATS,"MDS1")
plot_map(Test_STATS,"MDS2")
# Showing observed heterozygosity on map:
plot_map(Test_STATS,"Ho")
# Showing expected heterozygosity on map:
plot_map(Test_STATS,"Hs")
# Showing Fis on map:
plot_map(Test_STATS,"Fis")
