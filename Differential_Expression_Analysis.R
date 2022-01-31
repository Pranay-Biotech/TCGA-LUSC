
#Installing Bioconductor package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.12') #Check the recent version
#BiocManager::available() #To see abvailable packages in BioCManager
#BiocManager::install("GenomicFeatures") #To install any package out of displayed packages

library(TCGAbiolinks)
#Check GDC server status using the api https://api.gdc.cancer.gov/status
isServeOK() 

#To download gene Expression data
#Harmonized data (data aligned against genome of reference hg38)
query.exp <- GDCquery(project = "TCGA-LUSC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      experimental.strategy = "RNA-Seq",
                      workflow.type = "HTSeq - Counts",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))

#Legacy data ((data aligned against genome of reference hg19)
query.exp <- GDCquery(project = "TCGA-LUSC",
                      data.category = "Gene expression",
                      legacy=TRUE,
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq",
                      experimental.strategy = "RNA-Seq",
                      file.type ="normalized_results",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))

#For a given TCGA project it gets the samples (barcode) with both DNA methylation and Gene expression data from GDC database
#matchedMetExp("TCGA-LUSC",legacy=TRUE)

GDCdownload(query.exp,directory ="GDCdata",method="api",files.per.chunk=20)# can use method = "api", files.per.chunk = 10 (no. of chunks you want to divide your download)

View(as.data.frame(getResults(query.exp)))


#preparing the data 
exp<-GDCprepare(query.exp,save=TRUE,save.filename ="exp_prepared",directory="GDCdata")

library(SummarizedExperiment)
View(assay(exp)) #To see counts data
LuscMatrix<-assay(exp) #To get the matrix of gene counts

View(as.data.frame(colData(exp)))#To see metadata of samples 
sample_metadata<-(as.data.frame(colData(exp)))

View(as.data.frame((rowData(exp)))) #To see gene info
#gene_info<-as.data.frame(exp@rowRanges)

#Create a Summary table for each sample in a project saying if it contains or not files for a certain data category
sample_summary<-getDataCategorySummary("TCGA-LUSC", legacy = TRUE)
#View(getSampleFilesSummary("TCGA-LUSC"))


#Get hg19 or hg38 information from biomaRt (in case, it is not available in rowRanges df)
#get.GRCh.bioMart(genome = c("hg19", "hg38"), as.granges = FALSE)

#To visualise outliers
View(TCGAanalyze_Preprocessing(exp))


#If Counts dataframe has ensemble id's as identifiers, we have to replace it with gene codes
#Hence merging LiscMatrix with gene_info (in case of hg)



BiocManager::install("limma")
install.packages("limma")
browseVignettes("limma")
library(limma)


BiocManager::install("edgeR")

library(EdgeR)
library(limma)

library(DESeq2)

#TCGAanalyze_DEA allows user to perform Differentially expression analysis (DEA), using edgeR package or limma to identify differentially expressed genes (DEGs)
#dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(exp, geneInfo)
#dataFilt <- TCGAanalyze_Filtering(tabDF = dataBRCA, method = "quantile", qnt.cut =  0.25)
#samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                           # mat2 = dataFilt[,samplesTP],
                            #Cond1type = "Normal",
                            #Cond2type = "Tumor")
#TCGAanalyze_Normalization allows user to normalize mRNA transcripts and miRNA, using EDASeq package.
#TCGAanalyze_Normalization performs normalization using: 

#EDASeq::newSeqExpressionSet

#EDASeq::withinLaneNormalization

#EDASeq::betweenLaneNormalization

#EDASeq::counts


BiocManager::install("EDASeq")
library(EDASeq)
dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(exp, geneInfo)
View(as.data.frame(dataNorm))

#TCGAanalyze_Filtering allows user to filter mRNA transcripts and miRNA, samples, higher than the threshold defined quantile mean across all samples.
dataFilt <- TCGAanalyze_Filtering(tabDF =dataNorm , method = "quantile", qnt.cut =  0.25)

#Following are the arguments of TCGAanalyze_Filtering()
# tabDF	:
# is a dataframe or numeric matrix, each row represents a gene, each column represents a sample come from TCGAPrepare
#
# method	:
# is method of filtering such as 'quantile', 'varFilter', 'filter1', 'filter2'
# 
# qnt.cut	:
# is threshold selected as mean for filtering
# 
# var.func	:
# is function used as the per-feature filtering statistic. See genefilter documentation
# 
# var.cutoff:	
# is a numeric value. See genefilter documentation
# 
# eta	:
# is a parameter for filter1. default eta = 0.05.
# 
# foldChange :	
# is a parameter for filter2. default foldChange = 1.




#TCGAquery_SampleTypes return the union of samples that are from these type.
#TCGAquery_SampleTypes(barcode, typesample)
#NT :	Solid Tissue Normal
#TP	: PRIMARY SOLID TUMOR
samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

#TCGAtumor_purity Filters TCGA barcodes according to purity parameters
dataTableSubt <- TCGAtumor_purity(colnames(dataFilt),
                                  estimate = 0.6,
                                  absolute = 0.6,
                                  ihc = 0.8,
                                  lump = 0.8,
                                  cpe = 0.7)
#List with $pure_barcodes attribute as a vector of pure samples and $filtered attribute as filtered samples with no purity info
#57 pure barcodes and 51 filtered attributes 
samplesTP_pure <- TCGAquery_SampleTypes(dataTableSubt$pure_barcodes, typesample = c("TP"))
#we can also feed these pure_barcodes as tumor samples to TCGAanalyze_DEA fxn instead of samplesTP


#Differentially expression analysis (DEA) using edgeR or limma
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                             mat2 = dataFilt[,samplesTP],
                             Cond1type = "Normal",
                             Cond2type = "Tumor")
# there are Cond1 type Normal in  51 samples
# there are Cond2 type Tumor in  502 samples
# there are  6711 features as genes 


#Using pure_barcodes or samplesTP_pure as tumor samples obtained after selecting pure tumors
dataDEGs_pure<- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP_pure],
                            Cond1type = "Normal",
                            Cond2type = "Tumor")
#there are Cond1 type Normal in  51 samples
# there are Cond2 type Tumor in  57 samples
# there are  6711 features as miRNA or genes 



#dataDEG_limma<-TCGAanalyze_DEA_Affy(AffySet, FC.cut = 0.01)
#AffySet	: A matrix-like data object containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples



#Volcano Plot
TCGAVisualize_volcano(
  dataDEGs$logFC,
  dataDEGs$PValue,
  filename = NULL,
  y.cut = 10000000,
  x.cut = c(1),
  title= "Volcano plot (Normal vs Tumor)",
  legend = "Status",
  names.fill = TRUE,
  show = "both"
)
