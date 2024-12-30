# RADseq

RADseq data analysis for RIVERSCAPE project

### Read VCFs and sample data
### Calculate and save population parameters 
### Calculate and save Fst matrix for each species

```
library(readxl)
library(dplyr)
library(vegan)
library(rbiom)
library(tidyr)
library(vcfR)
library(ggplot2)
library(adegenet)
library(poppr)
library(dplyr)
library(hierfstat)
library(vegan)
library(leaflet)
library(mapmixture)
library(pheatmap)

 
#
setwd("F:/Riparian/ANALYSES_CLEAN/new_outfiles/")
setwd("/media/bartek/DANE/Riparian/ANALYSES_CLEAN/new_outfiles")

library(riverdist)
library(sp)

RiverNetwork <- line2network(path="F:/Riparian/ANALYSES_CLEAN/Carper", layer="RasterT_RIVERS_1_1",tolerance = 10,reproject="+init=epsg:25832")


all_results<-list()
VCF_NAMES<-list.files(path ="VCFs/")

FST_LIST<-list()
POPULATION_SUMMARY_LIST<-list()
loaded_genind$hierarchy

#  Iterating all VCFs to extract population summary and Fst, and
for(vcf_name in VCF_NAMES){

VCF_PATH<-paste0("VCFs/",vcf_name)
loaded_vcf<-read.vcfR(VCF_PATH) # loading vcf, for tests read part of the file using nrow=500
loaded_genind<-vcfR2genind(loaded_vcf) #conversion to gendind
data.hierfstat<-genind2hierfstat(loaded_genind,pop=gsub("-.*","",  gsub(".*_", "", rownames(loaded_genind$tab) ) ))
STATS<-basic.stats(data.hierfstat)
Ho_populations<-apply(STATS$Ho,2,function(x){mean(x,na.rm=T)})
Hs_populations<-apply(STATS$Hs,2,function(x){mean(x,na.rm=T)})
Fis_populations<-apply(STATS$Fis,2,function(x){mean(x,na.rm=T)})
Fst<-genet.dist(data.hierfstat) #,method="WC84")  #

METADATA<-read.csv2("SAMPLE_DATA.csv",sep=",")
METADATA<-METADATA[ METADATA$Sample %in% rownames(loaded_genind$tab),]
METADATA$Population<-gsub("^.*_(.*)", "\\1",  gsub("^(.*?)-.*", "\\1",METADATA$Sample))  
POPULATION_METADATA<-METADATA %>% group_by(Population) %>% summarise(Population=Population[1],Latitude=mean(as.numeric(Latitude)),Longitude =mean(as.numeric(Longitude)),N=n() )%>%data.frame    

POPULATION_METADATA$Ho<-Ho_populations[POPULATION_METADATA$Population]
POPULATION_METADATA$Hs<-Hs_populations[POPULATION_METADATA$Population]
POPULATION_METADATA$Fis<-Fis_populations[POPULATION_METADATA$Population]

write.table(POPULATION_METADATA,file=paste0(vcf_name,".population_summary.csv"))
write.table(as.matrix(Fst),file=paste0(vcf_name,".fst_matrix.csv"))


# EUCLIDEAN DISTANCES
euc_dist<-dist.mat(POPULATION_METADATA,"Population","Latitude","Longitude",  POPULATION_METADATA,"Population","Latitude","Longitude")
euc_dist<-data.frame(euc_dist)%>%dplyr::select(-from_to)
reshaped_euc_dist<-reshape(euc_dist, idvar = "from", timevar = "to", direction = "wide")
reshaped_euc_dist<-data.frame(reshaped_euc_dist)%>%dplyr::select(-from)
colnames(reshaped_euc_dist)<-POPULATION_METADATA$Population
rownames(reshaped_euc_dist)<-POPULATION_METADATA$Population

for(i in 1:ncol(reshaped_euc_dist)){
  reshaped_euc_dist[i,i]<-0
}
write.table(as.matrix(reshaped_euc_dist),file=paste0(vcf_name,".euclidean_distance_matrix.csv"))

# RIVER DISTANCES
SPATIAL_POINTS<-SpatialPointsDataFrame( coords = data.frame(Xbng =  POPULATION_METADATA$Longitude, Ybng =  POPULATION_METADATA$Latitude) ,data=POPULATION_METADATA ) 
proj4string(SPATIAL_POINTS) <- CRS("+init=epsg:4326") # WGS 84
SPATIAL_POINTS.transformed<- spTransform(SPATIAL_POINTS, "+init=epsg:25832") 
SPATIAL_POINTS.transformed.xy2segvert <- xy2segvert(  x= coordinates(SPATIAL_POINTS.transformed)[,1], y=coordinates(SPATIAL_POINTS.transformed)[,2], rivers= RiverNetwork)
hist(SPATIAL_POINTS.transformed.xy2segvert$snapdist, main="snapping distance (m)")

#riverpoints(seg=SPATIAL_POINTS.transformed.xy2segvert$seg, vert=SPATIAL_POINTS.transformed.xy2segvert$vert, rivers=RiverNetwork, pch=15,    col="blue")
#riverdistance(startseg=SPATIAL_POINTS.transformed.xy2segvert$seg[7], startvert=SPATIAL_POINTS.transformed.xy2segvert$vert[7],
#              endseg=SPATIAL_POINTS.transformed.xy2segvert$seg[15], endvert=SPATIAL_POINTS.transformed.xy2segvert$vert[15], rivers=RiverNetwork, map=TRUE)/1000

RIVER_DISTANCES<-riverdistancematbysurvey(indiv=1,unique=1, survey=SPATIAL_POINTS$Population, seg=SPATIAL_POINTS.transformed.xy2segvert$seg, 
                                                 vert=SPATIAL_POINTS.transformed.xy2segvert$vert, rivers=RiverNetwork,full=T)

write.table(as.matrix(RIVER_DISTANCES),file=paste0(vcf_name,".river_distance_matrix.csv"))
}

```


Plotting:
```
library(sf)
library(prettymapr)
library(ggspatial)
library(ggspatial)
library(rlang)
plot_map_ggplot <- function(data, var,title="") {
  # Ensure var is a column in the data
  if (!var %in% names(data)) {
    stop("Variable not found in data")
  }
  
  # Create a color scale based on the variable
  pal <- colorRampPalette(c("red", "green"))(100)
  color_scale <- scale_color_gradientn(colors = pal, limits = range(data[[var]], na.rm = TRUE))
 
  sf_data <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)
  

  # Create the map
  ggplot() +annotation_map_tile(type = "osm",zoom=12)+
     geom_sf(data = sf_data, aes(color = !!sym(var)), size = 4)  +
    color_scale +
    xlab("")+ylab("")+
    ggtitle(paste("Map of", var)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )+
    theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) +
      ggtitle(title)
}



PLOTS<-list()
for (j in 1:length(VCF_NAMES)){
XD<-read.table(paste0(VCF_NAMES[j],".population_summary.csv"),sep=' ')
p<-plot_map_ggplot(XD%>%filter(N>=5),"Hs", VCF_NAMES[j] )
PLOTS<-append(PLOTS,list(p))
}
PLOT_HS_ALL<-ggpubr::ggarrange(plotlist = PLOTS, ncol = 5, nrow = ceiling(length(PLOTS) / 5))
ggsave(PLOT_HS_ALL,file="HS_all_vcfs_new.pdf",width=20,height=20)


PLOTS<-list()
for (j in 1:length(VCF_NAMES)){
XD<-read.table(paste0(VCF_NAMES[j],".population_summary.csv"),sep=' ')
p<-plot_map_ggplot(XD,"Ho", VCF_NAMES[j] )
PLOTS<-append(PLOTS,list(p))
}
PLOT_HO_ALL<-ggpubr::ggarrange(plotlist = PLOTS, ncol = 5, nrow = ceiling(length(PLOTS) / 5))
ggsave(PLOT_HO_ALL,file="Ho_all_vcfs.pdf",width=20,height=20)


```



