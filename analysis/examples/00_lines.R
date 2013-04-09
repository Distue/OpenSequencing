# ------------------------------------------------------------
# Initialization of data
# ------------------------------------------------------------

source("/share/massstorage/seq-scripts/seqhelp.R")
project   <- "mycn-lines"
projectDir   <- "mycn-lines-mRNA-pe"
baseDir      <- file.path("/share/massstorage/projects/", projectDir)
dataDir      <- file.path(baseDir, "alignment/run2/")
phenoFile    <- file.path(baseDir, "raw/pheno.txt")
fastqDir     <- file.path(baseDir, "fastq/")
cufflinksDir <- file.path(dataDir, "cufflinks202") 
aligner      <- "tophat"
outputDir    <- file.path(dataDir, "routput")
pheno        <- read.delim(phenoFile, row.names=1)

orig.names.col <- "OrigName"
groups.col     <- "Cellline"

setwd(dataDir)
dir.create(outputDir, showWarnings=F)

myorder <- rownames(pheno)
myorder.short <- myorder

# create a list having the original sample names as values and the new names as list names
# the order defines the output order in the gene expression table
input.names            <- createInputNameList(pheno)

p.value = 0.01
# ------------------------------------------------------------
# Get Annotations from BioMART
# ------------------------------------------------------------

mappingLookupTable           <- getEnsemblMappingTable()
rownames(mappingLookupTable) <- mappingLookupTable[,1]

# ------------------------------------------------------------
# Create FPKM expression table
# ------------------------------------------------------------

# read in the cufflink files
#cufflinksFileList      <- readCufflinksFiles(input.names, cufflinksDir = cufflinksDir)
# create the expression table
# expressionTable        <- createExpressionTableFromCufflinks(cufflinksFileList    = cufflinksFileList,
#                                                             mappingLookupTable   = mappingLookupTable,
#                                                             mappingTableIdColumn = 0,
#                                                             action               = "sum")

# add sum and average FPKM column
#expressionTable        <- addStatsColumnsToExpressionTable(expressionTable, names(input.names))

# Check 
# expressionTable[duplicated(expressionTable[,"ids"]),1:2]

# write it as a text file
#write.table(as.matrix(expressionTable), 
#            file = file.path(cufflinksDir, paste(project, "ensembl.fpkm.gene.expression.table.txt", sep=".")),
#            sep="\t", quote=F, row.names=F)

# ------------------------------------------------------------
# Read FPKM expression table
# ------------------------------------------------------------

#expressionTable <- read.delim(file.path(cufflinksDir, paste(project, "ensembl.fpkm.gene.expression.table.txt", sep=".")), row.names=1, sep="\t", stringsAsFactors=F)

# ------------------------------------------------------------
# Summarization htseq
# ------------------------------------------------------------

htseqFiles <- file.path(dataDir, paste(paste(getColumnForIndex(orig.names.col, pheno), aligner, "sorted.sam.htseq.count.txt", sep=".")))
counts.unfiltered <- countTableFromHtseqCount(htseqFiles, column.names = rownames(pheno))

save(counts.unfiltered, file=file.path(outputDir, paste(project, "counts.Rdata", sep=".")))
write.table(counts.unfiltered, file=file.path(outputDir, paste(project, "counts.txt", sep=".")), sep="\t", quote=F)

# --------------------------------------------------------------------------------------
# Load the prestored counts
# --------------------------------------------------------------------------------------
load(file.path(outputDir, paste(project, "counts.Rdata", sep=".")))

# --------------------------------------------------------------------------------------
# getCPMKbs
# --------------------------------------------------------------------------------------
genelength <- getGeneLengthForCounts("/share/massstorage/seq-genomes/genelength.txt", counts=counts.unfiltered)
cpmkbs <- getCPMKBs(counts.unfiltered[,myorder], genelength)

cpmkbs.annotated <- annotateTable(as.data.frame(cpmkbs), 0, mappingLookupTable, 0, T)

write.table(cpmkbs.annotated, sep="\t", file=file.path(outputDir, paste(project, "CPMkb.txt", sep=".")), quote=F) 

cpmkbs.annotated <- read.table(sep="\t", file=file.path(outputDir, paste(project, "CPMkb.txt", sep="."))) 

# ------------------------------------------------------------
# Filter
# ------------------------------------------------------------

# filters counts for min of 1 CPM per minimum number of samples in smallest group
counts <- filterCounts(counts.unfiltered, getColumnForIndex(groups.col, pheno))

# ------------------------------------------------------------
# MDS
# ------------------------------------------------------------
require(edgeR)

plotMDS(counts)
plotMDS(counts[,c(1:4,7,8)])       
plotMDS(counts[,c(1:8)])       
 
# ------------------------------------------------------------
# Correlation
# ------------------------------------------------------------

cor(counts)
corcounts <- cor(counts)

combinations <- lapply(levels(pheno[,groups.col]), function(x) {
   combn(rownames(pheno[pheno[,groups.col] == x,]), m=2)[,1]
})
names(combinations) <- levels(pheno[,groups.col])
corlist <- lapply(combinations, function(x) { corcounts[x[[1]], x[[2]]]})

corcounts["KCN_A", "KCN_B"]
# 0.9988506
corcounts["KCNR_A", "KCNR_B"]
# 0.9975851
corcounts["Kelly_A", "Kelly_B" ]
# 0.9933007
corcounts["SY5Y_A", "SY5Y_B"]
# 0.9992131
corcounts["IMR32_A", "IMR32_B"]
# 0.9996092

# ------------------------------------------------------------
# Differential gene expression analysis 
# ------------------------------------------------------------

selection        <- c("SY5Y", "KCN", "KCNR", "Kelly", "IMR32")
pheno.selected   <- pheno[pheno[,groups.col] %in% selection,]
counts.selected  <- counts[,rownames(pheno.selected)]
cellline        <- factor(pheno.selected[,groups.col], levels = selection)

# edgeR with GLM function - tagwise dispersion.
dd       <- DGEList(counts = counts.selected, group = cellline)

design   <- model.matrix(~0 + cellline)  # Create a design matrix for calculating the GLM.
fit.tgw  <- edgeR.GLMfit(dd, design)

lrt <- list()
lrt[["SY5Y.vs.KCN"]]     <- glmLRT(fit.tgw, contrast =  c(-1,1,0,0,0))
lrt[["SY5Y.vs.KCNR"]]    <- glmLRT(fit.tgw, contrast =  c(-1,0,1,0,0))
lrt[["SY5Y.vs.Kelly"]]   <- glmLRT(fit.tgw, contrast =  c(-1,0,0,1,0))
lrt[["SY5Y.vs.IMR32"]]   <- glmLRT(fit.tgw, contrast =  c(-1,0,0,0,1))
lrt[["KCN.vs.KCNR"]]     <- glmLRT(fit.tgw, contrast =  c(0,-1,1,0,0))
lrt[["KCN.vs.Kelly"]]    <- glmLRT(fit.tgw, contrast =  c(0,-1,0,1,0))
lrt[["KCN.vs.IMR32"]]    <- glmLRT(fit.tgw, contrast =  c(0,-1,0,0,1))
lrt[["KCNR.vs.Kelly"]]   <- glmLRT(fit.tgw, contrast =  c(0,0,-1,1,0))
lrt[["KCNR.vs.IMR32"]]   <- glmLRT(fit.tgw, contrast =  c(0,0,-1,0,1))
lrt[["Kelly.vs.IMR32"]]  <- glmLRT(fit.tgw, contrast =  c(0,0,0,-1,1))

saveAllSignificantGenesAnnotate(lrt, cpmkbs.annotated, 0, project, outputDir, T, p.value=p.value)
writeLines(con=file.path(folder, paste(project, "-info.txt", sep="")), paste("p.value=", p.value, sep="")) 

# ------------------------------------------------------------
# Statistics: Amount of genes
# ------------------------------------------------------------


stats <- createOutputStatistic(list.files(outputDir, ".*lines.+vs.+txt$", full=FALSE),
                               outputDir)


write.table(stats, file=file.path(outputDir, paste(project, "stats.txt", sep=".")), sep="\t", quote=F)
            
# ------------------------------------------------------------
# DEseq filter for clustering
# ------------------------------------------------------------

# create special pheno file for original input from deseq
phenoDEseq <- pheno
phenoDEseq[,"filename"] <- paste(pheno[,orig.names.col], aligner, "sorted.sam.htseq.count.txt", sep=".") 
phenoDEseq[,"samplename"] <- rownames(pheno)
phenoDEseq <- phenoDEseq[,c("samplename", "filename", groups.col)]

cds <- newCountDataSetFromHTSeqCount(phenoDEseq, dataDir)                        
cds <- estimateSizeFactors(cds)
sizeFactors(cds)                         

# model formula has to be same as used in fits
cds <- estimateDispersions(cds, method="pooled-CR", modelFormula = count ~ condition)
vsd <- varianceStabilizingTransformation(cds)
p <- plotPCA(vsd, intgroup=c("condition"))

# specify full and complete models (fits)
fit1 <- fitNbinomGLMs(cds, count~condition)
fit0 <- fitNbinomGLMs(cds, count~1)

# test
pval <- nbinomGLMTest(fit1, fit0)
padj <- p.adjust(pval, method="BH")
res <- cbind(fit1, pval = pval, padj = padj)
head(res[order(res$padj),])

table(res$padj < 0.05)
table(res$padj < 0.01)


# filtered genes
resfiltered <- na.omit(res[res$padj < p.value,])
resfiltered <- resfiltered[rownames(resfiltered) %in% rownames(counts),]

# write those genes
write.csv(res, file=file.path(outputDir, paste(project, "res-DESeq.csv", sep="-")))
write.csv(resfiltered, file=file.path(outputDir, paste(project, "-res-DESeq-filtered-", p.value, ".csv", sep="")))
hist(res$pval, breaks=100)

# read the genes 
resfiltered <- read.csv(file=file.path(outputDir, paste(project, "-res-DESeq-filtered-", p.value, ".csv", sep="")), row.names=1)
head(resfiltered)                     

# get the adjusted expression values from the differentially expressed genes 
decluster <- na.omit(cpmkbs[rownames(resfiltered),myorder])
head(decluster)

# scale the genes
declusterscaled <- na.omit(t(scale(t(decluster))))


# ------------------------------------
# mclust
# ------------------------------------
library(mclust)
mclustfit <- Mclust(declusterscaled)
print(fit) # display the best model 
plot(fit) # plot results

# ------------------------------------
# kmeans
# ------------------------------------

tryKMeans(declusterscaled, 40)

fit <- kmeans(declusterscaled, 20, iter.max=100) # X cluster solution
# get cluster means
myclusters <- aggregate(declusterscaled, by=list(fit$cluster), FUN=mean)

printClusters(decluster, myorder, fit$cluster, declusterscaled, c(4,5), file.path(outputDir, "clustering-kmeans"))

# ------------------------------------
# hclust with correlation 
# ------------------------------------

# calculating the pearson correlation
all.genes.cor <- calculateFeatureCorrelations(decluster)

all.genes.cor[1:10,1:10]

# any NAs in there?
sum(is.na(all.genes.cor[1:10,]))

apply(all.genes.cor, 1, sd)

# create a distance matrix by subtracting pearson correlation coefficient from 1
all.genes.cor.dist <- as.dist(1 - all.genes.cor)

# do the hclust with average linkage
all.genes.tree <- hclust(all.genes.cor.dist, method='average')

plot(all.genes.tree, main=paste("Gene clustering ;", "Pearson distance", "; average linkage"), xlab=NULL, cex=0.1, cex.main=1.5)

# cut the tree
cluster.selection.hclust <- cutree(all.genes.tree, k=25)

# print clusters
printClusters(decluster, myorder, cluster.selection.hclust, declusterscaled, c(5,5), file.path(outputDir, "clustering-hclust"))

# overview how many genes overall in the clusters
hcluststat <- as.data.frame(table(cluster.selection))



