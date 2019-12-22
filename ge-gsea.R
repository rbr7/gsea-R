
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# foo <- listDatasets(mart)


f <- listFilters(mart)
# View(f)
# getBM(attributes = c("entrezgene", "description"), filters="hgnc_symbol", values = "IL10", mart=mart)
# or
# d <- getBM(attributes = c("entrezgene", "description", "refseq_mrna"), filters="refseq_mrna", values=E$genes$REFSEQ, mart=mart)

a <- listAttributes(mart)
# View(a)



# Matching existing REFSEQ ids to mrna refseq
d <- getBM(attributes=c("entrezgene_id", "description", "refseq_mrna"), filters="refseq_mrna", mart=mart, values=E$genes$REFSEQ)
En$genes$EntrezID <- d$entrezgene[ match(En$genes$REFSEQ, d$refseq_mrna) ]


# Repeat fit to get the EntrezIDs in fit object
d <- model.matrix(~ 0 + group, data=En$targets)
colnames(d) <- levels(En$targets$group)
c <- makeContrasts(TBvsNID="TB-NID", LTBIvsNID="LTBI-NID", TBvsLTBI="TB-LTBI", levels=d)
fit2 <- lmFit(En, d)
fit2 <- eBayes(contrasts.fit(fit2, c))

# running basic gsea
res <- goana(fit2, geneid=En$genes$EntrezID, coef="TBvsNID")
restable <- topGO(res, ontology="BP")
# View(restable)

# KEGG pathways analysis
res <- kegga(fit2, geneid=En$genes$EntrezID, coef="TBvsNID")
topKEGG(res)



# use topGO
# source("https://bioconductor.org/biocLite.R")
# biocLite("topGO")
library(topGO)

tt <- topTable(fit2, coef="TBvsNID", number=Inf, genelist=En$genes)  
tt <- tt[ !is.na(tt$EntrezID), ]   
tt <- tt[!duplicated(tt$EntrezID), ]  
rownames(tt) <- tt$EntrezID  


universe <- tt$adj.P.Val   
names(universe) <- tt$EntrezID
data("geneList")

# caculating topology
go <- new("topGOdata", ontology="BP", 
  allGenes=universe,
  geneSel=topDiffGenes, nodeSize=10, 
  annotationFun=annFUN.org, 
  mapping="org.Hs.eg.db")

# tests
resF    <- runTest(go, algorithm="classic", statistic="fisher" ) 
resKS   <- runTest(go, algorithm="classic", statistic="ks" )   
resF.e  <- runTest(go, algorithm="elim", statistic="fisher" )  
resKS.e <- runTest(go, algorithm="elim", statistic="ks" ) 
GenTable(go, classicFisher=resF, classicKS=resKS, elimFisher=resF.e,
             elimKS=resKS.e)

# using GOrilla
write.csv(tt$GENE_SYMBOL, row.names=FALSE, quote=FALSE, file="export.csv")
# go to http://cbl-gorilla.cs.technion.ac.il/ and upload the file
# filter based on 


# tmod
tt <- topTable(fit2, coef="TBvsNID", number=Inf, sort.by="p")
gnames <- tt$GENE_SYMBOL

# install.packages("tmod")
library(tmod)
res <- tmodCERNOtest(gnames)
head(res)


# View(res)
# View(tt)
# showModule(tt$REFSEQ, x$GENE_SYMBOL, module = "LI.M7.0") 
# foo <- showModule(tt, x$GENE_SYMBOL, module = "LI.M7.0") 
# View(foo)

data(tmod)
names(tmod)
# View(tmod$MODULES)
tmod$MODULES2GENES$LI.M7.0


# inspect res
head(res)

# inspect single modules
evidencePlot(tt$GENE_SYMBOL, m = res$ID[1:3], col=2:4, legend = "right")  
# evidencePlot(gnames, m = res$ID[1:3], col=2:4, legend = "right")  
evidencePlot(gnames, m = "LI.M7.0", col=2:4, legend = "right", gene.labels = TRUE) 

# CD96 
showGene(En$E, En$targets$group) 
tt[ tt$GENE_SYMBOL] 

showGene(En$E[ "A_23_P44155", ], En$targets$group)

library(tagcloud)
tagcloud(res$Title, weights=-log10(res$P.Value), col=smoothPalette(res$AUC))
## size of tags=P value, color=AUC


res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL)
# res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=go.bp)
tmodPanelPlot(res, text.cex=0.6, grid="between", filter.rows.pval = 1e-9)


# tmod : panel plot
tmodPanelPlot(res, text.cex=0.7, filter.rows.pval = 1e-7)


pie <- tmodLimmaDecideTests(fit2, genes=fit2$genes$GENE_SYMBOL, 
                            pval.thr = 0.01, lfc.thr = 0.5)
tmodPanelPlot(res, text.cex=1.1, filter.rows.pval = 1e-7, pie=pie)


tmodPanelPlot(res, text.cex=0.6, filter.rows.pval = 1e-7, pie=pie,
              pie.style = "r", grid = "between")

# Using MSigDB database

# Go to: http://software.broadinstitute.org/gsea/downloads.jsp 
# Select the most recent data for "All gene sets" under current MSigDB xml file
msig <- tmodImportMSigDB("msigdb_v6.1.xml")

# testing for modules
# E.g. KEGG pathways
kegg <- msig[ msig$MODULES$Subcategory == "CP:REACTOME" ]  
# kegg <- msig[ msig$MODULES$Subcategory == "CP:KEGG" ]  


res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=kegg)

hallmark <- msig[ msig$MODULES$Category == "H" ]
hallmark$MODULES$Title <- gsub("Hallmark", "", hallmark$MODULES$Title)
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=hallmark)
# C7 
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=msig$MODULES$Category == "C7") 


with(msig$MODULES, table(Category, Subcategory))

go.bp <- msig[ msig$MODULES$Subcategory == "BP" ]
res <- tmodLimmaTest(fit2, fit2$genes$GENE_SYMBOL, mset=go.bp)

# Using KEGG
library(KEGGREST) 
pathways <- as.matrix(keggLink("pathway", "hsa"))
pathways <- keggLink("pathway", "hsa")
paths    <- sapply(unique(pathways), function(p) keggGet(p)[[1]]$NAME) 
paths <- gsub(" - Homo sapiens.*", "", paths)

# building tmod dataset 
m2g <- sapply(unique(pathways), function(p) names(pathways)[pathways == p ], simplify=F)
head(m2g)
m <- data.frame(ID=unique(pathways), Title=paths)
kegg2 <- makeTmod(modules=m, modules2genes=m2g)
res <- tmodLimmaTest(fit2, paste0("hsa:", fit2$genes$EntrezID), mset=kegg2)

# KEGGREST for pathway
png <- keggGet("path:hsa05322", option="image")
grid::grid.raster(png)

# source("https://bioconductor.org/biocLite.R")
# biocLite("pathview")
library(pathview)
install.packages("png")
library(png)
genes <- tt$logFC
names(genes) <- tt$EntrezID
foo <- pathview(genes, species="hsa", pathway.id="05322")
image <- readPNG("hsa05322.pathview.png")
grid::grid.raster(image)



sizes <- seq(3, 25, by=3)
tbs <- which(En$targets$group == "TB")
nids <- which(En$targets$group == "NID")
sel <- sapply(sizes, function(s) {
  c(sample(tbs, s), sample(nids, s))
}, simplify=FALSE)
names(sel) <- paste0("S.", sizes) 

 
# for various sample sizes
res <- sapply(sel, function(s) {
  ee <- En[ , s]
  ee$targets$group <- factor(ee$targets$group)
  d <- model.matrix(~ group, data=ee$targets)
  f <- eBayes(lmFit(ee, d))
  tt.temp <- topTable(f, number=Inf, p.value=0.05)
  print("Number of significant genes:")
  print(nrow(tt.temp))
  res <- tmodLimmaTest(f, f$genes$GENE_SYMBOL, coef=2)
  res <- res[[1]]
}, simplify=FALSE)

tmodPanelPlot(res, filter.rows.pval =1e-7)
