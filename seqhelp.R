
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
                                               idColumn="ids",
                                               trackingidColumn="tracking_id",
                                               locusColumn="locus",
                                               FPKMColumn="FPKM") { 

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
   
   if(!is.list(cufflinksFileList)) {
      stop("cufflinksFileList is not a list")
   }
   
   # create a list which contains all identifiers
   all.identifiers <- c()
   
   
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
   expression.table            <- data.frame(row.names=all.identifiers)
   expression.table[,idColumn] <- unlist(lapply(all.identifiers, function(x) return(strsplit(x, "\\|")[[1]][[1]])))
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
   
   # get rid of duplicates by adding | 1 or | 2 at the end of the id
   expression.table <- handleDuplicateRows(expression.table, idColumn=idColumn, action=action, samplenames=names(cufflinksFileList))
   
   # remove the id column
   # rownames(expression.table) <- expression.table[,idColumn]
   # expression.table <- expression.table[,-which(colnames(expression.table) == idColumn)]
   
   return(expression.table)
}




# ------------------------------------------------------------
# Summarization
# ------------------------------------------------------------
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

annotateTable <- function(dataTable, colTable, mappingTable, colMapping, expression.matrix = NULL) {
   # vector which stores the ids for annotation
   idsForAnnotation <- c()
   # if colTable is 0, rownames are taken
   if (colTable == 0) {
      idsForAnnotation <- rownames(dataTable)
   } else { # else the column stored in colTable is taken  
      idsForAnnotation <- dataTable[,colTable]
   }
   
   resultTable <- t(sapply(idsForAnnotation, function(x) {
      # select rows to annotate
      annotationRows <- as.data.frame(mappingTable[mappingTable[,colMapping] == x,])
      
      # check if there are more than 2 annotations
      if(nrow(annotationRows) > 1) {
         # warning("Warning, multiple Mapping problem")
         tmp <- t(apply(annotationRows, 2, function(y) {
            as.character(paste(unique(y), collapse=", "))
         }))
         
         rownames(tmp) <- x
         
         return(tmp)
      } else {
         return(annotationRows[1,])
      }
   }))
   
   if (is.null(expression.matrix)) {
      return(cbind(resultTable, dataTable))
   }
   else {
      lapply(resultTable[,2], function(x) {
         sum(rownames(expression.matrix) == x)
         expression.matrix[x,]
      })
      
      return(cbind(resultTable, dataTable))
   }
}


# extracts the table of significant genes from GLM
getSignificantGenes <- function(lrtObject, adjust.method = "fdr") {
   # get the count of the up and downregulated genes and sum them up
   countDGE <- sum(summary(decideTestsDGE(lrtObject, adjust.method=adjust.method))[c("-1","1"),])
   
   topTagsObj <- topTags(lrtObject, n = countDGE, adjust.method=adjust.method)
   
   topTagsTable <- topTagsObj$table 
   colnames(topTagsTable)[ncol(topTagsTable)] <- topTagsObj$adjust.method
   
   return(topTagsTable)
}


getSignificantGenesAndAnnotate <- function(lrtObject, mappingTable, expression.matrix = NULL) { 
   dataTable <- getSignificantGenes(lrtObject)  
   return(annotateTable(dataTable = dataTable, colTable = 0, mappingTable = mappingTable, expression.matrix = expression.matrix, colMapping = 1))
}


saveSignificantGenesAndAnnotate <- function(lrtObject, mappingTable, expression.matrix = NULL, name, folder) {
   tmptable <- getSignificantGenesAndAnnotate(lrtObject, mappingTable = mappingTable, expression.matrix = expression.matrix)
   write.table(as.matrix(tmptable), file=file.path(folder, paste(name, ".txt", sep="")), sep="\t", quote=F, row.names=F)
}


# ------------------------------------------------------------
# Normalize Counts
# ------------------------------------------------------------

normalizeCounts <- function(counts, pheno) {
   cond <- as.factor(pheno[colnames(counts),2])
   cds <- newCountDataSet(counts, cond)
   cds <- estimateSizeFactors( cds )  
   return(counts(cds, normalized=TRUE))
}
