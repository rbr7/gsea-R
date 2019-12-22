
# reference article
# Data from Tuch, B.B. et al. (2010). Tumor transcriptome sequencing
# reveals allelic expression imbalances associated with copy number 
# alterations. PLoS ONE 5, e9317
rawdata <- read.csv("rnaseq_example.csv", stringsAsFactors=FALSE)

library(edgeR)
library(statmod)
y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])  


# filtering data
keep <- rowSums(cpm(y) > 0.5) > 1
dim(y)
y <- y[keep, ]
dim(y)


# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")

library(org.Hs.eg.db)
eg2refseq <- toTable(org.Hs.egREFSEQ)
y <- y[ y$genes$idRefSeq %in% eg2refseq$accession, ] 
y$genes$EntrezID <- eg2refseq$gene_id[ match(y$genes$idRefSeq, eg2refseq$accession) ] 

# head(y$genes)

# removing duplicates
head(which(duplicated(y$genes$EntrezID)))
ord <- order(rowSums(y$counts), decreasing=TRUE)  
y <- y[ord, ]
y <- y[ !duplicated(y$genes$EntrezID), ]

# calculating normalization factors
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
plotMDS(y)



x <- y$counts
vars <- apply(y$counts, 1, var)
keep <- vars > 0 
x <- x[keep,]
x <- cpm(x, log=T, prior.count=3)

# pca
pca <- prcomp(t(x), scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19, col=rep(c(1,2), 3))


# preparing model
Patient <- factor(paste0("P.", c(8,8,33,33,51,51)))
Tissue <- factor(c("N","T","N","T","N","T"))
design <- model.matrix(~Patient+Tissue)
# View(design)

y <- estimateDisp(y, design, robust=TRUE)

# fitting model
fit <- glmQLFit(y, design)

# calculate results for : coefficient of interest
lrt <- glmQLFTest(fit, coef="TissueT")
topTags(lrt)

# significant genes
de <- decideTestsDGE(lrt, p.value=0.01, lfc=1)
table(de)

#install GO.db
#biocLite("GO.db")
library(GO.db)

go <- goana(lrt)
topGO(go, ont="BP", sort="Up", n=30)


tt <- topTags(lrt, n=Inf)
write.csv(tt$table$nameOfGene, file="output.txt", row.names = FALSE, quote = FALSE)

# Use the GOrilla web site -> upload "output.txt"
