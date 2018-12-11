library(Matrix)
library(pagoda2)
library(dplyr)
markers <- read.csv("/Data/celltype_genesets.txt", header = F, sep = '\t', na.strings=c("","NA"), stringsAsFactors = F)
cm <- readRDS(file.path(find.package('pagoda2'),'extdata','sample_BM1.rds'))
dim(cm)
cm[1:3,1:3]
par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0)
hist(log10(colSums(cm)+1),main='molecules per cell',col='cornsilk',xlab='log10(molecules per cell)')
hist(log10(rowSums(cm)+1),main='molecules per gene',col='cornsilk',xlab='log10(molecules per gene])')
counts <- gene.vs.molecule.cell.filter(cm,min.cell.size=500)
hist(log10(rowSums(counts)+1),main='Molecules per gene',xlab='molecules (log10)',col='cornsilk')
abline(v=1,lty=2,col=2)
counts <- counts[rowSums(counts)>=10,]
dim(counts)
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts,log.scale=TRUE, n.cores=2)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=50,n.odgenes=3e3)
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=igraph::infomap.community,type='PCA')
M <- 30; r$getEmbedding(type='PCA',embeddingType = 'largeVis', M=M,perplexity=30,gamma=1/M,alpha=1)
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,min.group.size=50,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (largeVis)')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=F)
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (tSNE)')
r$getKnnClusters(method=igraph::multilevel.community,type='PCA',name='multilevel')
r$getKnnClusters(method=igraph::walktrap.community,type='PCA',name='walktrap')
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',groups=r$clusters$PCA$community,show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='infomap clusters (tSNE)')
r$plotEmbedding(type='PCA',embeddingType='tSNE',clusterType='multilevel',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='multlevel clusters (tSNE)')
r$getDifferentialGenes(type='PCA',verbose=T,clusterType='community')
de <- r$diffgenes$PCA[[1]][['2']];
r$plotGeneHeatmap(genes=rownames(de)[1:15],groups=r$clusters$PCA[[1]])
suppressMessages(library(org.Hs.eg.db))
str(de)
str(r$counts)
r$counts[1:3, 1:10]
ids <- unlist(lapply(mget(colnames(r$counts),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
str(rids)
rids[1:10]
rids <- names(ids); names(rids) <- ids;
go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
hdea <- r$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='community',z.threshold=3)
str(hdea)
hdea[1:2]
gene_sets <- hierDiffToGenesets(hdea)
str(head(genesets))
library(GO.db)
termDescriptions <- Term(GOTERM[names(go.env)]); # saves a good minute or so compared to individual lookups
head(names(go.env))
sn <- function(x) { names(x) <- x; x}
gene_sets.go <- lapply(sn(names(go.env)),function(x) {
  list(properties=list(locked=T,genesetname=x,shortdescription=as.character(termDescriptions[x])),genes=c(go.env[[x]]))
})
gene_sets <- c(gene_sets, gene_sets.go)
deSets <- get.de.geneset(r, groups = r$clusters$PCA[['community']], prefix = 'de_')
gene_sets.custom <- apply(markers,1,function(x) {
  list(properties=list(locked=T,genesetname=as.character(x[1]),shortdescription=as.character("No description")),genes=unique(as.vector(x[c(-1,-2)][as.vector(x[c(-1,-2)][!is.na(x[(c(-1,-2))])]) %in% cm@Dimnames[[1]]])))
})
str(head(gene_sets.custom))
gene_sets <- c(gene_sets, gene_sets.custom)
appmetadata <- list(apptitle = 'October_Demo_App')
r$makeGeneKnnGraph(n.cores = 2)
additionalMetadata <- list()
additionalMetadata$community <- p2.metadata.from.factor(r$clusters$PCA[['community']], displayname = 'Infomap', s = 0.7, v = 0.8,start = 0.1, end = 0.5)
additionalMetadata$multilevel <- p2.metadata.from.factor(r$clusters$PCA[['multilevel']], displayname = 'Multilevel', s = 0.9, v = 0.8,start = 0.5, end = 1)
a <- r$clusters$PCA[['walktrap']]
library(colorRamps)
p1 <- colorRamps::primary.colors(n = nlevels(a))
names(p1) <- levels(a)
additionalMetadata$walktrap <- p2.metadata.from.factor(r$clusters$PCA[['walktrap']], displayname = 'Walktrap', pal = p1)
p2web <-
  make.p2.app(
    r,
    dendrogramCellGroups = r$clusters$PCA$community,
    additionalMetadata = additionalMetadata,
    geneSets = gene_sets.custom,
    appmetadata = appmetadata,
    show.clusters = FALSE # Hide the clusters that were used for the dendrogram from the metadata
  );
show.app(app=p2web,name='app')
