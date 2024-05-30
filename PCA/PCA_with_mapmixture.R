# Plotting PCA from VCF using mapmixture

# input files:
vcf<-"input_files/test.vcf"
sample_data<-read.csv("input_files/test_sample_data.csv")

# output files:
output_filename<-"output_files/PCA_1.pdf"




library(mapmixture)
library(ggplot2)
library(adegenet)
library(RColorBrewer)
library(gridExtra)
library(vcfR)
 
loaded_vcf<-read.vcfR(vcf)
loaded_genind<-vcfR2genind(loaded_vcf)
pop(loaded_genind) <-as.vector(sample_data$Population)
loaded_genind<-filtered_genind

# Converting to hierfstat to deal with missing data by indpca function
Data.hierfstat<-genind2hierfstat(loaded_genind)
PCA<-indpca(Data.hierfstat)

# calculation percent of explained variance by axes
percent = round(PCA$ipca$eig/sum(PCA$ipca$eig)*100, digits = 2)

scatter_pca <- scatter_plot(
  dataframe = PCA$ipca$li,
  group_ids = loaded_genind$pop,
  type = "points",
  axes = c(1,2),
  percent = percent,
  #colours = cols,
  point_size = 3,
  point_type = 21,
  centroid_size = 4,
  stroke = 0.1,
  plot_title = "PCA coloured by populations"
)+
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 10),
  )
# plotting
scatter_pca
# saving
ggsave(output_filename,scatter_pca)




