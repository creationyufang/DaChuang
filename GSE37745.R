##Step 1 Getting Started

#install the required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("packages you want")
library(GEOquery)
library(Biobase)
library(multiClust)
library(preprocessCore)
library(ctc)
library(gplots)
library(dendextend)
library(graphics)
library(grDevices)
library(amap)
library(survival)

gse <- getGEO(filename="GSE37745_series_matrix.txt.gz")
data.gse <- exprs(gse)
pheno <- pData(phenoData(gse))
WriteMatrixToFile(tmpMatrix=data.gse, tmpFileName="GSE37745.expression.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)
WriteMatrixToFile(tmpMatrix=pheno, tmpFileName="GSE37745.clinical.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)

#Normalization of Gene Expression Datasets
data.norm <- normalize.quantiles(data.gse, copy=FALSE)

# shift data before log scaling to prevent errors from log scaling negative numbers
if (min(data.norm)> 0) {
  print("perfect data")
} else {
  mindata.norm=abs(min(data.norm)) + .001
  data.norm=data.norm + mindata.norm  
}

# Log2 scaling of the dataset
data.log <- t(apply(data.norm, 1, log2))

# Write the gene expression and clinical data to text files
WriteMatrixToFile(tmpMatrix=data.log,
                  tmpFileName="GSE37745.normalized.expression.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)


#$$Formatting the Patient Clinical Information
#before using system.file function, please mv the files(expression and clinical) to
# system.file("extdata", "GSE2034.normalized.expression.txt", package= "multiClust")
clin_file <- system.file("extdata", "GSE37745-RFS-clinical-outcome.txt",
                         package="multiClust")
clinical <- read.delim2(file=clin_file, header=TRUE)
clinical[1:5, 1:2]




##Step 2 Loading Your Gene Probe Expression Dataset into R

exp_file <- system.file("extdata", "GSE37745.normalized.expression.txt", package= "multiClust")
data.exprs <- input_file(input=exp_file)
data.exprs[1:4,1:4]


##Step 3 Gene Selection Algorithms

#Choosing 50% of the total selected gene probes in a dataset
gene_num <- number_probes(input=exp_file, 
                          data.exp=data.exprs, Fixed=NULL,
                          Percent=50, Poly=NULL, Adaptive=NULL)
#method = CV_Rank  CV_Guided SD_Rank Poly
ranked.exprs <- probe_ranking(input=exp_file,
                              probe_number=gene_num, 
                              probe_num_selection="Fixed_Probe_Num",
                              data.exp=data.exprs, 
                              method="SD_Rank")


##Step 4 Cluster Analysis of Selected Genes and Samples

#The gap_statistic option has a very long computational time 
# and can take up to several hours!
# so I choose Fixed Cluster Number
cluster_num <- number_clusters(data.exp=data.exprs, Fixed=3,
                               gap_statistic=NULL)

# Call the cluster_analysis function
#cluster_type=Kmeans or HClust
kmeans_analysis <- cluster_analysis(sel.exp=ranked.exprs,
                                    cluster_type="Kmeans",
                                    distance=NULL, linkage_type=NULL, 
                                    gene_distance=NULL, num_clusters=cluster_num,
                                    data_name="GSE37745", probe_rank="SD_Rank",
                                    probe_num_selection="Fixed_Probe_Num",
                                    cluster_num_selection="Fixed_Clust_Num")


##Step 5 Obtaining the Average Expression for Each Gene/Probe in Each Cluster

avg_matrix <- avg_probe_exp(sel.exp=ranked.exprs,
                            samp_cluster=kmeans_analysis,
                            data_name="GSE37745", cluster_type="Kmeans", distance=NULL,
                            linkage_type=NULL, probe_rank="SD_Rank",
                            probe_num_selection="Fixed_Probe_Num",
                            cluster_num_selection="Fixed_Clust_Num")
head(avg_matrix)



##Step 6 Clinical Analysis of Selected Gene Probes and Samples

# Call the avg_probe_exp function
surv <- surv_analysis(samp_cluster=kmeans_analysis, clinical=clin_file,
                      survival_type="RFS", data_name="GSE37745", 
                      cluster_type="Kmeans", distance=NULL,
                      linkage_type=NULL, probe_rank="SD_Rank",
                      probe_num_selection="Fixed_Probe_Num", 
                      cluster_num_selection="Fixed_Clust_Num")
