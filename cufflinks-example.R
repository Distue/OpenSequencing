# Created by: Thomas Schwarzl <thomas@schwarzl.net>
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# script for creating an expression table from cufflinks
# --------------------------------------------------------------------------------------
source("/share/massstorage/seq-scripts/seqhelp.R")

# --------------------------------------------------------------------------------------
# input 
# --------------------------------------------------------------------------------------

dataDir      <- "/share/massstorage/projects/mycn-sy5y-mRNA-pe/alignment/run1/"
cufflinksDir <- file.path(dataDir, "cufflinks202") 
pheno        <- read.delim("/share/massstorage/projects/mycn-sy5y-mRNA-pe/raw/pheno.txt", row.names=1)
setwd(dataDir)

# --------------------------------------------------------------------------------------
# create FPKM expression table 
# --------------------------------------------------------------------------------------

# create mapping table with one IDs as row names (mappingTableIdColumn = 0) which is fast
# or as column (mappingTableIdColumn = columnname/columnID) which is slower
mappingLookupTable           <- getEnsemblMappingTable()
rownames(mappingLookupTable) <- mappingLookupTable[,1]

# create a list having the original sample names as values and the new names as list names
# the order defines the output order in the gene expression table
input.names            <- createInputNameList(pheno)

# read in the cufflink files
cufflinksFileList      <- readCufflinksFiles(input.names, cufflinksDir = cufflinksDir)

# create the expression table
expressionTable        <- createExpressionTableFromCufflinks(cufflinksFileList    = cufflinksFileList,
                                                             mappingLookupTable   = mappingLookupTable,
                                                             mappingTableIdColumn = 0,
                                                             action               = "sum")

# write it as a text file
write.table(as.matrix(expressionTable), 
            file = file.path(cufflinksDir, "ensembl.fpkm.gene.expression.table.txt"),
            sep="\t", quote=F, row.names=F)

