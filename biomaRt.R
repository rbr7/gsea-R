
# load packages
library(biomaRt)
View(listMarts())


mart <- useMart("ensembl")


# loading datasets for organism
View(listDatasets(mart))

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


View(listAttributes(mart))
View(listFilters(mart))


bm <- getBM(filters="illumina_humanwg_6_v3", 
            values = kaforou$genes$ID, 
            attributes = c( "hgnc_symbol", "description", "reactome", "illumina_humanwg_6_v3", "go_id", "name_1006", "go_linkage_type"), 
            mart=mart)

# probe ids
sum(kaforou$genes$ID %in% bm$illumina_humanwg_6_v3)
all(bm$illumina_humanwg_6_v3 %in% kaforou$genes$ID)

# unique go ids
goids <- unique(bm$go_id)
length(goids)

# mapping
go2gene <- split(bm$illumina_humanwg_6_v3, f = bm$go_id)

head(lengths(go2gene))
max(lengths(go2gene))

go2gene <- go2gene[ lengths(go2gene) >= 20 ]
go2gene <- go2gene[ lengths(go2gene) <= 250 ]
length(go2gene)

goterms <- data.frame(ID=names(go2gene), 
   Title=bm$name_1006[ match(names(go2gene), bm$go_id)])

# reverse mapping : gene2go
gene2go <- split(rep(names(go2gene), lengths(go2gene)), unlist(go2gene))
mygo <- makeTmod(modules=goterms, modules2genes = go2gene, genes2modules = gene2go)

res <- tmodLimmaTest(fit, genes=fit$genes$ID, mset=mygo)
