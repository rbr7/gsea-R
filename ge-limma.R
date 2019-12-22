
# loading limma
library(limma)


# E <- new("EListRaw", list(E=E.mat, genes=genes, targets=ph))


plotDensities(E, log=TRUE, legend=FALSE)


# E.bg <- backgroundCorrect(E, method="normexp")
# plotDensities(E.bg, log=TRUE, legend=FALSE)


# normalizing between arrays
En <- normalizeBetweenArrays(E, method="quantile")
plotDensities(En, log=TRUE, legend=FALSE)


# rm(E, E.bg, geo)


dim(En)


# unique probes
length(unique(En$genes$NAME))

# averaging
En <- avereps(En, ID=En$genes$NAME) 
dim(En)


En <- En[ En$genes$CONTROL_TYPE == "FALSE", ]
dim(En)


# pca 
pca <- prcomp(t(En$E), scale.=TRUE)  
plot(pca$x[,1], pca$x[,2], pch=19)


# pca 3d view
library(pca3d)
pca3d(pca)

# pca for genes
pca <- prcomp(En$E, scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19)
pca3d(pca)


pca <- prcomp(t(En$E), scale.=TRUE)
plot(pca$x[,1], pca$x[,2], pch=19, col=factor(En$targets$group))
legend("topright", as.character(levels(factor(En$targets$group))), pch=19, col=1:3)
plot(pca$x[,3], pca$x[,4], pch=19, col=factor(En$targets$group))
legend("topright", as.character(levels(factor(En$targets$group))), pch=19, col=1:3)


# important genes
ord <- order(abs(pca$rotation[,3]), decreasing = T)
head(En$genes[ord,])


# filtering
iqrs <- apply(En$E, 1, IQR)
cutoff <- quantile(iqrs, 0.90) 
sum(iqrs < cutoff) / nrow(En) * 100 
En.f <- En[ iqrs > cutoff, ] 



En$targets$group <- factor(En$targets$group, levels=c("NID", "LTBI", "TB"))
d <- model.matrix(~ group, data=En$targets)
# View(d) 

fit1 <- lmFit(En, d)
fit1 <- eBayes(fit1)
tt <- topTable(fit1, coef="groupTB") 
# View(tt)


# effect of filtering
#fit1b <- eBayes(lmFit(En.f, d))

#nrow(topTable(fit1, coef="groupTB", number=Inf, p.value=0.001, lfc=1))
#nrow(topTable(fit1b, coef="groupTB", number=Inf, p.value=0.001, lfc=1))

#tt1 <- topTable(fit1, coef="groupTB", number=Inf, sort.by="none")
#tt1b <- topTable(fit1b, coef="groupTB", number=Inf, sort.by="none")
#tt1 <- tt1[ rownames(tt1b), ] # only genes present in tt1b as well
#plot(tt1$logFC, tt1b$logFC)
#plot(tt1$P.Value, tt1b$P.Value, log="xy", pch=19, col="#33333311")
#abline(0,1, col="red")
#filtering did not have a impact


d <- model.matrix(~ 0 + group, data=En$targets)
# View(d)

colnames(d) <- levels(En$targets$group)
fit2 <- lmFit(En, d)

c <- makeContrasts(TBvsNID="TB-NID", LTBIvsNID="LTBI-NID", TBvsLTBI="TB-LTBI", 
                   TBvsAll="TB-(NID+LTBI)/2", 
                   levels=d)
fit2 <- contrasts.fit(fit2, c)
fit2 <- eBayes(fit2) # bayes factors, p-value
# head(fit1$coefficients) 
# head(fit2$coefficients)  

# results
cor(fit1$coefficients[,3], fit2$coefficients[,1])
smoothScatter(fit1$coefficients[,3], fit2$coefficients[,1]) 


topTable(fit2, coef="TBvsNID")



volcanoplot(fit2, coef="TBvsNID")


tt <- topTable(fit2, coef="TBvsNID", number=Inf)
with(tt, plot(logFC, -log10(P.Value), pch=19))
# top 50 genes
with(tt[1:50,], points(logFC, -log10(P.Value), pch=19, col="red"))
with(tt[1:50,], text(logFC, -log10(P.Value), labels=GENE_SYMBOL, col="red"))



En$targets$gr.sex <- paste0(En$targets$group, ".", En$targets$sex)
En$targets$gr.sex <- factor(En$targets$gr.sex)

d <- model.matrix(~ 0 + gr.sex, data=En$targets)
colnames(d) <- levels(En$targets$gr.sex)
fit3 <- lmFit(En, d)

c <- makeContrasts(int="(TB.Female-NID.Female)-(TB.Male-NID.Male)", levels=d)
fit3 <- eBayes(contrasts.fit(fit3, c))
topTable(fit3, coef=1)


# heatmaps
tt <- topTable(fit2, coef="TBvsNID", number=Inf)
sel <- rownames(tt)[1:40]

# select data
x <- data.matrix(En$E[sel,])
library(gplots)
genelabs <- En$genes$GENE_SYMBOL[ match(sel, En$genes$NAME)]
heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group)


colf <- colorRampPalette(c("purple", "black", "cyan"))
heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group, col=colf)

heatmap.2(x, trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group, col=colf, 
             breaks=seq(-2, 2, length.out=21))

# remove extra information
ord <- order(En$targets$group)
heatmap.2(x[,ord], trace="n", labRow = genelabs, scale="r", 
             labCol=En$targets$group[ord], col=colf, 
             breaks=seq(-2, 2, length.out=21),
             key = F, Colv = NULL, colsep=62)
