BiocManager::install("topGO")
BiocManager::install("Rgraphviz")

library(topGO)
library(Rgraphviz)

# Set wd
setwd("gene enrichment analysis/nut morphology traits")

# Genes associated to nut traits analysed all together

### Load Data
annotation <- readMappings(file = "CavTom2PMs-1.0_gene_universe.txt")
geneUniverse <- names(annotation) 
genesOfInterest <- read.table("list_all_genes_nut_morphology_traits.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

### GO Enrichment Analyses

# BP
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = annotation)
myGOdata 
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata) 
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher") 
resultFisher
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)
allRes
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "list_all_genes_nut_morphology_traits_bp", pdfSW = TRUE)

# CC
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = annotation)
myGOdata 
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata) 
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher") 
resultFisher
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)
allRes
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "list_all_genes_nut_morphology_traits_cc", pdfSW = TRUE)

# MF
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = annotation)
myGOdata 
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata) 
resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher") 
resultFisher
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)
allRes
showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')
printGraph(myGOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "list_all_genes_nut_morphology_traits_mf", pdfSW = TRUE)

