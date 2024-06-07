# Plotting PCA from VCF using mapmixture
#
# 2 functions:
# PCA_from_vcf() -  Calculates PCA using path to VCF file, and path to CSV file with population assignment
# plot_PCA() -  Plots PCA based on an object returned by PCA_from_vcf function



library(mapmixture)
library(ggplot2)
library(adegenet)
library(RColorBrewer)
library(gridExtra)
library(vcfR)
library(dplyr)
library(hierfstat)

PCA_from_vcf<-function(vcf_path,sample_data_path){
  # Loading sample metadata
  sample_data<-read.csv(sample_data_path)
  # Loading vcf file
  loaded_vcf<-read.vcfR(vcf_path)
  #conversion to genind
  loaded_genind<-vcfR2genind(loaded_vcf)
  
  # joining sample info with sample names from VCF
  samples_joined_with_info<-left_join( data.frame(Sample=rownames(loaded_genind$tab)), sample_data)
  
  #adding info on population assignment to individuals
  pop(loaded_genind) <-as.vector(samples_joined_with_info$Population)
  
  ## Add filtering here:
  ## filtered_genind <- filter(loaded_genind) 
  ## loaded_genind<-filtered_genind
  ##
  
  
  # Converting to hierfstat to deal with missing data by indpca function
  Data.hierfstat<-genind2hierfstat(loaded_genind)
  PCA<-indpca(Data.hierfstat)
  
return(PCA)
}


# plot_PCA -  Plots PCA based on an object returned by PCA_from_vcf function
plot_PCA<-function(PCA,axes_X=1,axes_Y=2,title="PCA coloured by populations"){

# calculation percent of explained variance by axes
percent = round(PCA$ipca$eig/sum(PCA$ipca$eig)*100, digits = 2)

scatter_pca <- scatter_plot(
  dataframe = PCA$ipca$li,
  group_ids = factor(PCA$ipca$rownames),
  type = "points",
  axes = c(axes_X,axes_Y),
  percent = percent,
  #colours = cols,
  point_size = 3,
  point_type = 21,
  centroid_size = 4,
  stroke = 0.1,
  plot_title = title
)+
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 10),
  )

return(scatter_pca)

}



########################################################################################################
# Tests:

setwd("F:/Riparian/ANALYSES_CLEAN/PCA")
vcf_path_carper<-"../Carper/Carper.vcf"
sample_data_path_carper<-"../Carper/Carper_sample_data.csv"
pca_object_carper<-PCA_from_vcf(vcf_path_carper,sample_data_path_carper)
pca_plot_carper<-plot_PCA(pca_object_carper,axes_X = 1,axes_Y = 2,title = "Carduus personata")

vcf_path_alninc<-"../Alninc/Alninc.vcf"
sample_data_path_alninc<-"../Alninc/Alninc_sample_data.csv"
pca_object_alninc<-PCA_from_vcf(vcf_path_alninc,sample_data_path_alninc)
pca_plot_alninc<-plot_PCA(pca_object_alninc,axes_X = 1,axes_Y = 2,title = "Alnus incaca")

ggpubr::ggarrange(pca_plot_alninc,pca_plot_carper,ncol=2)

# saving
ggsave("output_files/Alninc_PCA_1.pdf",pca_plot_alninc)
ggsave("output_files/Carper_PCA_1.pdf",pca_plot_carper)


