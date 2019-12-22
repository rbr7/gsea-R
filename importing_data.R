
# new directory and change the workspace
library(this)
this()


# load packages
library(GEOquery)
geo <- getGEO(GEO="GSE28623") # getting GEO dataset will take time


geo
class(geo)
str(geo)
names(geo)


geo <- geo[[1]]


# Converting the "ExpressionSet" object into an "EList" object from limma
ph <- as(phenoData(geo), "data.frame")
E.mat  <- exprs(geo)
genes <- as(featureData(geo), "data.frame")


View(ph)


head(ph)
ph <- ph[ , c("geo_accession", "source_name_ch1", "characteristics_ch1.1")]
colnames(ph) <- c("ID", "group", "sex")
ph$group <- gsub(".* ", "", ph$group)  
ph$sex <- gsub(".* ", "", ph$sex) 
head(ph)


# Actual expression data
str(E.mat)
head(E.mat[,1:10])


colnames(genes)


# select genes
genes <- genes[ , c("NAME", "CONTROL_TYPE", "REFSEQ", "GENE_SYMBOL", "DESCRIPTION")]
head(genes)
tail(genes)


# fitting
dim(E.mat)
dim(genes)
dim(ph)


# loading limma to create an EList
library(limma)
E <- new("EListRaw", list(E=E.mat, genes=genes, targets=ph))
str(E)
class(E)
E[1:10,1:2]
E[1:10,1:4]

dim(E)

