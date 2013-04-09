require(WriteXLS)
require(edgeR)
require(DESeq)
require(gplots)
# -----------------------------------------------------------------
# Cufflinks to Expression Table
# -----------------------------------------------------------------


# create a list with the original names of the samples and the new names
# input:   pheno - a dataframe containing the pheno information
#          sampleNewName - defines column to be taken as NewNames
#                          0 means rownames, number is the column position, 
#                          character string is the column name
#          sampleOrigNames - defines column which contains the original names
#                            same as properties as sampleNewNames
# hint - if the names stay the same, just create a list yourself with the 
#        sample names as values and listnames
createInputNameList <- function(pheno, sampleNewNames=0, sampleOrigNames="OrigName") {
   input.names <- list()
   
   # new names as list
   if(sampleOrigNames == 0) {
      input.names <- as.list(rownames(pheno))
   } else {
      input.names <- as.list(as.character(pheno[,sampleOrigNames]))
   }
   
   # orig names as names
   if (sampleNewNames == 0) {
      names(input.names) <- rownames(pheno)
   } else {
      names(input.names) <- as.character(pheno[,sampleNewNames])
   }
   
   return(input.names)
}


# Read Cufflinks Files
# input: input.names - list containing the cufflinks output directories and the names of the list represent
#                      the names how the samples should be called in the expression table
#        cufflinksDir - directory where the cufflinks output dirs are located
#                     - default: current working directory
#        cufflinksFileNames - what files to be taken:
#                           - default: "genes.fpkm_tracking"

readCufflinksFiles <- function(input.names, cufflinksDir=".", cufflinksFileName = "genes.fpkm_tracking") {
   cufflinksFiles <- list()
   if (!is.list(input.names)) {
      warning("input.names is not a list")
   }
   
   for(sample in names(input.names)) {
      inputFileName <- file.path(cufflinksDir, input.names[[sample]], cufflinksFileName)
      if(!file.exists(inputFileName)) {
         warning(paste("file ", inputFileName, " does not exist", sep=""))
      } else { 
         cufflinksFiles[[sample]]   <- read.table(inputFileName, header=T)
      }
   }
   
   return(cufflinksFiles)
}


# Expression Table from Cufflinks File List
# input: cufflinksFileList - list generated from readCufflinksFiles
#        action - rewrite: rewrites the id of duplicated rows, adds |1 |2 .. |n to the end of the id
#                 remove: removes the duplicated rows
#                 average: averages the rows
#                 median: takes the median of the rows
#                 sum: sums the rows up
#                 keepmax: keeps the row with the maximum sum
#        mappingLookupTable - is optional, if provided, the ids in trackingidColumn will be mapped
#                 with the ids in mappingTableIdColumn and the content will be added to the expression
#                 table.
#        mappingTableIdColumn - 0 is rownames, 1,2, .. n specifies the column number, and a character
#                 string would be the column name
#        idColumn - how the column with the ids should be names
#        trackingidColumn - where to find the identifier in the cufflinks file
#        locusColumn - where to find the locus in the cufflinks file
#        FPKMColumn - where to find the FPKMs in the cufflinks file
createExpressionTableFromCufflinks <- function(cufflinksFileList, action=c("sum", "rewrite", "remove", "keepmax", "average"),
                                               mappingLookupTable=F, mappingTableIdColumn = 0,
                                               idColumn = "ids",
                                               trackingidColumn = "tracking_id",
                                               locusColumn = "locus",
                                               FPKMColumn = "FPKM") { 


   
   if(!is.list(cufflinksFileList)) {
      stop("cufflinksFileList is not a list")
   }
   
   # create a list which contains all identifiers
   all.identifiers <- c()
   
   # get all identifiers
   for(i in names(cufflinksFileList)) {   
      # since the gene name can be duplicated - multiple locations for example 
      # the locus is needed to make the location unique
      # therefore the rownames for the table are combined id of tracking_id and locus, separated by | 
      rownames(cufflinksFileList[[i]]) <- paste(as.character(cufflinksFileList[[i]][,trackingidColumn]),
                                                as.character(cufflinksFileList[[i]][,locusColumn]), sep="|")   
      
      # store the (almost) unique identifier [TRACKINGID|LOCUS] in all.identifiers
      all.identifiers <- c(all.identifiers, rownames(cufflinksFileList[[i]]))
   }
   
   # only get unique identifiers 
   all.identifiers <- unique(sort(all.identifiers))
   
   # create expression matrix using the unique identifiers
   expression.table                <- data.frame(row.names=all.identifiers)
   expression.table[,idColumn]     <- unlist(lapply(all.identifiers, function(x) return(strsplit(x, "\\|")[[1]][[1]])))
   expression.table[,locusColumn]  <- unlist(lapply(all.identifiers, function(x) return(strsplit(x, "\\|")[[1]][[2]])))
   
   # if there is mapping information available, add it
   if(!is.logical(mappingLookupTable))  {
      # append the annotation
      expression.table <- cbind(expression.table, t(sapply(expression.table[,idColumn], function(x) {
         if(mappingTableIdColumn == 0) {
            return(mappingLookupTable[x,])
         }
         else {
            mappingLookupTable[mappingLookupTable[, mappingTableIdColumn] == x,]
         }
      })))
   }
   
   # fill in the values for the expression matrix
   # this is a very slow process
   for(i in names(cufflinksFileList)) {
      #expression.table[, i] <- 0
      expression.table[all.identifiers, i] <- cufflinksFileList[[i]][all.identifiers,FPKMColumn]
   }
   
   # all NAs -> 0
   expression.table[is.na(expression.table)] <- 0
   
   # get rid of duplicates 
   expression.table <- handleDuplicateRows(expression.table, idColumn=idColumn, action=action, samplenames=names(cufflinksFileList))
   
   # remove the id column
   # rownames(expression.table) <- expression.table[,idColumn]
   # expression.table <- expression.table[,-which(colnames(expression.table) == idColumn)]
   
   return(expression.table)
}


getDuplicatedRowNumbers <- function(fulltable, id, idColumn) {
   return(which(fulltable[,idColumn] == id))
}

rewriteDuplicateRow  <- function(fulltable, id, idColumn) {
   # get the rownumbers for the duplicated rows
   rownumbers <- getDuplicatedRowNumbers(fulltable, id, idColumn)
   
   # get the rows
   i <- 1
   for(j in rownumbers) {
      fulltable[ j, idColumn ] <- paste(  fulltable[ j, idColumn ], i, sep="|")
      i <- i + 1
   }
   
   return(fulltable)
}


keepmaxDuplicateRow <- function(fulltable, id, idColumn, samplenames) {
   rownumbers <- getDuplicatedRowNumbers(fulltable, id, idColumn)
   
   rows <- fulltable[rownumbers,]
   # get the sum of the rows
   sum.of.rows <- apply(rows[,samplenames], 1, sum)
   # decide which to keep
   row.to.keep <- which(sum.of.rows == max(sum.of.rows))
   
   if(length(row.to.keep) != 1) {
      stop(paste("could not determine the maximum of rows for ", id, sep=" "))
   }
   
   return(rbind(fulltable[-rownumbers,], rows[row.to.keep,]))
}

changeDuplicateRow <- function(fulltable, id, idColumn, samplenames, method) {
   rownumbers <- getDuplicatedRowNumbers(fulltable, id, idColumn)
   
   rows <- fulltable[rownumbers,]
   # calculate the new row based on method
   new.row.data <- apply(rows[,samplenames], 2, method)
   # restores metadata
   new.row <- rows[1,!colnames(rows) %in% samplenames]
   new.row[,names(new.row.data)] <- new.row.data
   rownames(new.row) <- paste(new.row[,idColumn], paste("changed", "to", deparse(substitute(method)), sep="_"), sep="|")
   return(rbind(fulltable[-rownumbers,], new.row))
}

# handles duplicated rows
handleDuplicateRows  <- function(fulltable, idColumn, action, samplenames) {
   # get duplicated rows   
   duplications <- as.character(sort(fulltable[duplicated(fulltable[,idColumn]),idColumn]))
   
   # for all duplications, do action
   switch(action, 
          rewrite = {
             for(id in duplications) {
                fulltable <- rewriteDuplicateRow(fulltable, id, idColumn=idColumn)
             }
          },
          remove = {
             fulltable <- fulltable[!fulltable[,idColumn] %in% duplications,]
          },
          average = {
             for(id in duplications) {
                fulltable <- changeDuplicateRow(fulltable, id, idColumn=idColumn, samplenames=samplenames, method=mean)
             }            
          },
          median = {
             for(id in duplications) {
                fulltable <- changeDuplicateRow(fulltable, id, idColumn=idColumn, samplenames=samplenames, method=median)
             }            
          },
          sum = {
             for(id in duplications) {
                fulltable <- changeDuplicateRow(fulltable, id, idColumn=idColumn, samplenames=samplenames, method=sum)
             } 
          },
          keepmax = {
             for(id in duplications) {
                fulltable <- keepmaxDuplicateRow(fulltable, id, idColumn=idColumn, samplenames=samplenames)
             }           
          },
         {
            stop(paste("action ", action, " not valid in handleDuplicateRows", sep=""))
         }
          
   )
   
   return(fulltable)
}

addStatsColumnsToExpressionTable <- function(expressionTable, sample.names) {
   expressionTable[,"FPKMsum"] <- apply(expressionTable[,sample.names], 1, sum)
   expressionTable[,"FPKMmean"] <- apply(expressionTable[,sample.names], 1, mean)
   expressionTable[,"FPKMsd"] <- apply(expressionTable[,sample.names], 1, sd)
   return(expressionTable)
}

# ------------------------------------------------------------
# Count Table from HTseq Count
# ------------------------------------------------------------

countTableFromHtseqCount <- function(htseqFiles, column.names) {
   htseqTables <- list()
   for(i in htseqFiles) {
      htseqTables[[i]] <- read.table(i, row.names=1)  
   }
   
   counts <- htseqTables[[1]]
   
   for( i in 2:length(htseqTables)) {
      counts <- cbind(counts, htseqTables[[i]])
   }
   
   colnames(counts) <- column.names
   
   counts <- counts[!rownames(counts) %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"),]
}


# -----------------------------------------------------------
# Filtering
# -----------------------------------------------------------

filterCounts <- function (counts, groups) {
   # calculate smallest group size
   smallestGroupSize <- min(summary(groups))
   
   # caluclate counts per million
   cpms   <- cpm(counts)
   # keep only the rows which more than 1 transcript than smallest group
   keep   <- rowSums(cpms > 1) >= smallestGroupSize
   
   # return selected counts
   counts <- counts[keep,]
}


# ------------------------------------------------------------
# Summarization
# ------------------------------------------------------------

# Usage:
# ------------------------------------------------------------
# bamFls <- list.files(dataDir, "sorted.bam$", full=TRUE)
# counts <- doSummarization(bamFls)
# 
# 
# colnames(counts) <- unlist(lapply(colnames(counts), function(x) {
#    rownames(pheno[pheno[,orig.names.col] == x,])
# }))
# 
# writeLines(paste("Writing counts to", outputDir, "counts.rData.", sep=" ")) 
# save(counts, file=file.path(outputDir, "counts.rData"))


getSampleNames <- function(coln, pheno, orig.names.col) {
   return(unlist(lapply(colnames(counts), function(x) {
      rownames(pheno[pheno[,orig.names.col] == x,])
   })))
}
                           

# for human hg19
doSummarization <- function(bamFls) {
   dir.create(outputDir, showWarnings = T)
   require(GenomicRanges)
   require(ShortRead)
   
   # fls <- list.files(fastqDir, "fastq$", full=TRUE)
   # names(fls) <- sub(".fastq", "", basename(fls))
   # ## use FastqSampler if fastq files are large
   # qas <- lapply(seq_along(fls),
   #               function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
   #               fls)
   # qa <- do.call(rbind, qas)
   # save(qa, file=file.path(outputDir, "qa.rda")
   #      browseURL(report(qa))
   #      
   #      
   
   # install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene")
   require(TxDb.Hsapiens.UCSC.hg19.knownGene)
   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
   # txdb <- Hsapiens_UCSC_hg19_ensGene_TxDb
   gnModel <- exonsBy(txdb, "gene")
   
   
   names(bamFls) <- sub("\\..*", "", basename(bamFls))
   counter <- function(fl, gnModel) {
      aln <- readGappedAlignments(fl)
      strand(aln) <- "*" # for strand-blind sample prep protocol
      hits <- countOverlaps(aln, gnModel)
      counts <- countOverlaps(gnModel, aln[hits==1])
      names(counts) <- names(gnModel)
      counts
   }
   
   counts <- sapply(bamFls, counter, gnModel)
   
   return(counts)
}


edgeR.GLMfit <- function(dd, design) {
   dd       <- estimateGLMCommonDisp(dd, design)
   dd       <- estimateGLMTrendedDisp(dd, design)
   dd       <- estimateGLMTagwiseDisp(dd, design)
   return(glmFit(dd, design))
}

# ------------------------------------------------------------
# Functions to annotate gene lists 
# ------------------------------------------------------------

library(biomaRt)

# function for the generation of a mapping table using BioMART, default: human
createMappingTable <- function(attributes, filters, dataset = "hsapiens_gene_ensembl", values = TRUE, uniqueRows = TRUE) {
   # query biomaRt and return results 
   return(getBM(attributes = attributes,
                filters    = filters,
                values     = values,
                mart       = useMart("ensembl", dataset = dataset),
                uniqueRows = uniqueRows
   ))
}

getEnsemblMappingTable <- function() {
   return(createMappingTable(
      attributes = c("ensembl_gene_id", "external_gene_id", "description"),
      filters = ""
   ))
}

getEntrezMappingTable <- function() {
   return(createMappingTable(
      attributes = c("entrezgene", "external_gene_id", "description"),
      filters = "with_entrezgene"
   ))
}

annotateTable <- function(dataTable, colTable, mappingTable, colMapping, appendBack=F, sep=", ") {
   checkIndex(colTable, dataTable)
   checkIndex(colMapping, mappingTable)
   
   # vector which stores the ids for annotation
   idsForAnnotation <- getColumnForIndex(colTable, dataTable)

   resultTable <- as.data.frame(t(vapply(idsForAnnotation, function(x) {
      # select rows to annotate
      annotationRows <- as.data.frame(getEntryForIndex(colMapping, mappingTable, x))

      return(as.character(t(apply(annotationRows, 2, function(y) {
         as.character(paste(unique(y), collapse=sep))
      })), stringsAsFactors=F))
   
   }, as.character(getEntryForIndex(colMapping, mappingTable, 1)) )), stringsAsFactors=F)
   
   colnames(resultTable) <- colnames(mappingTable)
   # append to the end or to the beginning
   if(appendBack == FALSE) {
      return(cbind(resultTable, dataTable))
   }
   else {
      return(cbind(dataTable, resultTable))
   }
}


# extracts the table of significant genes from GLM
getSignificantGenes <- function(lrtObject, adjust.method = "BH", p.value=0.05) {
   # get the count of the up and downregulated genes and sum them up
   countDGE <- sum(summary(decideTestsDGE(lrtObject, adjust.method=adjust.method, p.value=p.value))[c("-1","1"),])
   
   topTagsObj <- topTags(lrtObject, n = countDGE, adjust.method=adjust.method)
   
   topTagsTable <- topTagsObj$table 
   colnames(topTagsTable)[ncol(topTagsTable)] <- topTagsObj$adjust.method
   
   return(topTagsTable)
}


getSignificantGenesAndAnnotate <- function(lrtObject, mappingTable, colMapping, appendBack, p.value) { 
   checkIndex(colMapping, mappingTable)
   
   dataTable <- getSignificantGenes(lrtObject, p.value=p.value)  
   return(annotateTable(dataTable = dataTable, colTable = 0, mappingTable = mappingTable, colMapping = colMapping, appendBack = appendBack))
}

saveAllSignificantGenesAnnotate <- function(lrtList, mappingTable, colMapping, project, folder, appendBack, p.value=0.05, prefix=F) {
   if(!is.character(prefix)) {
      prefix <- ""
   } else {
      prefix <- paste("-", prefix, sep="")
   }
   
   require(WriteXLS)
   xlsxfile <- file.path(folder, paste(project, "-diff-genes.xlsx", sep=""))

   dfs <- list()
   for(i in names(lrtList)) {
      dfs[[i]] <- saveSignificantGenesAndAnnotate(lrtList[[i]], mappingTable, colMapping, name=paste(project, prefix, "-", i, sep=""), folder=folder, appendBack=appendBack, p.value=p.value)
   }

   WriteXLS("dfs", ExcelFileName=file.path(folder, paste(project, prefix, "-diff-genes.xls", sep="")),
               SheetNames=names(dfs), row.names=F, col.names=T, AutoFilter=T, BoldHeaderRow=T,
               FreezeRow = 1)
}

saveSignificantGenesAndAnnotate <- function(lrtObject, mappingTable, colMapping, name, folder, appendBack, p.value) {
   tmptable <- getSignificantGenesAndAnnotate(lrtObject, mappingTable = mappingTable, colMapping = colMapping, appendBack = appendBack, p.value = p.value)
   saveTable <- cbind(rownames(tmptable), tmptable)
   colnames(saveTable)[[1]] <- "ID"
   write.table(saveTable, file=file.path(folder, paste(name, ".txt", sep="")), sep="\t", quote=F, row.names=F) 
   return(saveTable)
}


checkIndex <- function(index, df, iname = deparse(substitute(index)), dfname = deparse(substitute(df))) {
   checkDataFrame(df, dfname=dfname)
   
   if(!(is.numeric(index) || is.character(index))) {
      stop(paste(index, "is not a number or character"))
   }  
   if(is.numeric(index)) {
      if(!(index >= 0 && index <= ncol(df))) {
         stop(paste(index, "is not 0 or within the numbers of columns of", dfname, sep=" "))
      }
   }
}

checkDataFrame <- function(df, dfname = deparse(substitute(df))) {
   if(!is.data.frame(df)) {
      stop(paste(dfname, "is not a data frame"))
   }   
}
   
getColumnForIndex <- function(i, df) {
   if (i == 0) { # if 0 rownames are taken
      return(rownames(df))
   } else { # else the column stored in colTable is taken  
      return(df[,i])
   }
}

getEntryForIndex <- function(i, df, x) {
   if (i == 0) { # if 0 rownames are taken
      return(df[x,])
   } else { # else the column stored in colTable is taken  
      return(df[df[,i] == x,])
   }
}

   
# ------------------------------------------------------------
# Normalize Counts
# ------------------------------------------------------------

normalizeCounts <- function(counts, cond) {
   cds <- newCountDataSet(counts, cond)
   cds <- estimateSizeFactors( cds )  
   return(counts(cds, normalized=TRUE))
}

# ------------------------------------------------------------
# Stats
# ------------------------------------------------------------

createOutputStatistic <- function(file.names, directory, logFCs=c(0.5,1,2,3)) {
   stats <- data.frame(row.names=file.names)
   
   for( i in file.names) {
      # read in the file
      tmptable <- read.delim(file.path(outputDir, i), header=T)
      
      # all
      stats[i,"all"] <- nrow(tmptable)

      for(fc in logFCs) {
         stats[i,paste(">",fc, "logFC", sep="")] <- nrow(tmptable[tmptable[,"logFC"] > fc,])
         stats[i,paste("<",fc, "logFC", sep="")] <- nrow(tmptable[tmptable[,"logFC"] < (-1 * fc),])
      }
   }
   
   return(stats)
}

# ------------------------------------------------------
# RPKM/FPKM calculation 
# ------------------------------------------------------

getGeneLengthForCounts <- function(file, counts) {
   genelength <- read.table(file, sep="\t", row.names=1)
   genelength <- genelength[rownames(counts),]
   return(genelength)
}

getUpperQuantileFPKMs <- function(counts, genelength) { 
   upperquantile <- apply(counts, 2, function(x) as.numeric(quantile(x, probs=0.75)))
   kUFPKMs <- counts / upperquantile / genelength * 1e+9 
   return(kUFPKMs)
}

getKiloUpperQuantileFPKMs <- function(counts, genelength) { 
   return(getUpperQuantileFPKMs(counts, genelength) / 1000)
}

getUpperQuantileFPKMs <- function(counts, genelength) { 
   upperquantile <- apply(counts, 2, function(x) as.numeric(quantile(x, probs=0.75)))
   UFPKMs <- counts / upperquantile / genelength * 1e+9
   return(UFPKMs)
}


getMegaMedianFPKMs <- function(counts, genelength) { 
   upperquantile <- apply(counts, 2, median)
   MMFPKMs <- counts / upperquantile / genelength * 1e+9 / 100000
   return(MMFPKMs)
}


getFPKMs <- function(counts, genelength) {
   totallibrary <- apply(counts, 2, sum)
   FPKMs <- counts / totallibrary / genelength * 1e+9   
}

getCPMBs <- function(counts, genelength) {
   return(cpm(counts) / genelength)
}

getCPMKBs <- function(counts, genelength) {
   return(cpm(counts) / (genelength/1e3))
}

getLog10CPMKBs <- function(counts, genelength, offset=1e3) {
   return(log10(getCPMKBs(counts, genelength) + offset) - log10(offset))
}


# ------------------------------------------------------
# Hierachical clustering
# ------------------------------------------------------

# to cope with the standard deviation is zero problem
filterForCor <- function(input) {
   non.zero.var <- logical()
   for(i in 1:ncol(input)) {
      non.zero.var[i] <- var(input[,i]) > 0
   }
   return(input[,non.zero.var])
}

calculateFeatureCorrelations <- function(data, use="pairwise.complete.obs", method="pearson") {
   return(cor(filterForCor(t(as.matrix(data))), use=use, method=method))
}

# ------------------------------------------------------
# K-means clustering
# ------------------------------------------------------

tryKMeans <- function(data, n.tries = 40, iter.max = 40, ...) {
   # find within group ss for all the data
   wss1 <- (nrow(data) - 1) * sum(apply(data, 2, var))
   wss  <- numeric(0)
   
   # calculate within group ss for 2:6 group partitions given by k-means clustering
   for(i in 2:n.tries) {
      W <- sum(kmeans(data, i, iter.max)$withinss)
      wss <- c (wss, W)
   }
   
   wss <- c(wss1, wss)
   plot (1:n.tries,
         wss, 
         type = "l", 
         xlab = "Number of groups", 
         ylab = "Within groups sum of squares", 
         lwd = 2,
         #xaxs = n.tries,
         ...)
}

# ------------------------------------------------------
# Clustering print Clusters
# ------------------------------------------------------

plotClusterOverview <- function(myclusters, cluster.selection, i) {
   cluster.genes <- names(cluster.selection[cluster.selection==i])
   plot(1:ncol(myclusters), myclusters[i,], type="o", ann=FALSE, axes=F)
   title(main=paste("Cluster-profile ", i, " (", length(cluster.genes), " genes)", sep=""))
   axis(1, at=1:ncol(myclusters), lab=sub("_", " ", colnames(myclusters)), las=2)
   
}

plotClusterGenes <- function(myclusters, cluster.selection, i, declusterscaled, myorder)  {
   cluster.genes <- names(cluster.selection[cluster.selection==i])
   plot(1:ncol(myclusters), myclusters[i,], type="o", ann=FALSE, axes=F, col="white", ylim=c(min(declusterscaled[cluster.genes,]),max(declusterscaled[cluster.genes,])))
   title(main=paste("Cluster ", i, " (", length(cluster.genes), " genes)", sep=""))
   axis(1, at=1:ncol(myclusters), lab=sub("_", " ", myorder), las=2)
   #axis(2, at=)
   for (j in na.omit(cluster.genes)) lines(1:ncol(myclusters), declusterscaled[j,myorder])
   
}

plotHeatmap <- function(declusterscaled, cluster.selection, i, myorder) {
   cluster.genes <- names(cluster.selection[cluster.selection==i])
   if(length(cluster.genes) > 1) {
      heatmap.2(as.matrix(declusterscaled[cluster.genes,myorder]),
                col     = colorRampPalette(c("darkblue", "white","darkred"))(100),
                scale   = "none",
                dendrogram = "row",
                trace   = "none", 
                hclust  = function(x) hclust(x,method="complete"),
                distfun = function(x) as.dist((1-cor(t(x)))/2),
                labCol  = sub("_", " ", myorder),
                main    = paste("Cluster ", i, " (", length(cluster.genes), " genes)", sep=""),
                key     = TRUE,
                keysize = 1,
                density.info = "density",
                Colv    = F
                
      )
   }
}

printClusters <- function(decluster, myorder, cluster.selection, declusterscaled, parWindows, outputDir, mappingTable) {
   dir.create(outputDir, showWarnings=F)
   myclusters <- aggregate(decluster[,myorder], by=list(cluster.selection), FUN=median)[,-1]

   # plot cluster overview 
   png(file.path(outputDir, paste("All-clusters-overview",".png", sep="")),  height=900, width=900, units="px")
   par(mfrow=parWindows) #, oma=c(0,0,0,0), mar=c(1,1,1,1))
   for(i in 1:nrow(myclusters)) {
      plotClusterOverview(myclusters, cluster.selection, i)
   }
   dev.off()
   
   # plot the clusters
   for(i in 1:nrow(myclusters)) {
      png(file.path(outputDir, paste("cluster", i, "overview.png", sep="-")))
      plotClusterOverview(myclusters, cluster.selection, i)
      dev.off()
      
      png(file.path(outputDir, paste("cluster", i, "allgenes.png", sep="-")))
      plotClusterGenes(myclusters, cluster.selection, i, declusterscaled, myorder)
      dev.off()
      
      #par(oma=c(1,1,1,1), mar=c(0,0,0,0))
      png(file.path(outputDir, paste("cluster", i, "heatmap.png", sep="-")))
      plotHeatmap(declusterscaled, cluster.selection, i, myorder)
      dev.off()
      
      df <- mappingTable[names(cluster.selection[cluster.selection==i]), ]
      WriteXLS("df", ExcelFileName=file.path(outputDir, paste("cluster", i, "allgenes.xls", sep="-")),
               SheetNames=paste("cluster", i, sep="-"), row.names=F, col.names=T, AutoFilter=T, BoldHeaderRow=T,
               FreezeRow = 1)
   }
}



write.xlstable <- function(x, file, names) {
   WriteXLS("x", ExcelFileName=file,
            SheetNames=names, row.names=F, col.names=T, AutoFilter=T, BoldHeaderRow=T,
            FreezeRow = 1)  
}

write.datatable <- function(x, file) {
   tmpTable <- cbind(rownames(x), x)
   colnames(tmpTable)[[1]] <- "ID"
   write.table(tmpTable, sep="\t", row.names=F, file=file, quote=F)
}
