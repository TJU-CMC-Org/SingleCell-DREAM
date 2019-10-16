rm(list = ls())
# Change the paths of the input files to rerun the script ###
normalized.data <- read.csv("Data/dge_normalized.txt.gz", sep = "\t") # 
# Change the paths of the input files to rerun the script ###
geometry <- read.csv("Data/geometry.txt.gz",sep = " ")

# Load datasets
load("Results_Common/my_mapCells_run.RData")

# Extract locations for RNAseq data
mccScores <- dm@mcc.scores
colnames(mccScores) <- colnames(normalized.data)
dim(mccScores)

# Index of cells position RNAseq position
index.RNAseq <- apply(X = mccScores, MARGIN = 2, FUN = function(coli){
    max_i <- max(coli)
    return(which(coli == max_i))
})

# There are MORE THAN ONE unique positions per cell
hist(sapply(X = index.RNAseq, FUN = length), 
     xlab = "Number of cell positions", 
     main = "287 Cells are NOT uniquely mapped", las = 1)

# For the cells with more than one positions take the first one
index.RNAseq.first <- sapply(X = index.RNAseq, 
                             FUN = function(li){li[[1]]})
geometry.RNAseq.all <- geometry[index.RNAseq.first, ]
geometry.RNAseq.all$Cell <- colnames(normalized.data)

# Keep the cells which are uniquely mapped
index.RNAseq.unique <- unlist(index.RNAseq[sapply(X = index.RNAseq, FUN = length) == 1])
geometry.RNAseq.unique <- geometry[index.RNAseq.unique, ]
geometry.RNAseq.unique$Cell <- names(index.RNAseq.unique)

# Save the true positions
save(geometry.RNAseq.all, geometry.RNAseq.unique, file = "Results_Common/geometry.RNAseq.RData")

