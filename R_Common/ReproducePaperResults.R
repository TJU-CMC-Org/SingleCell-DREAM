# Rerunning distMap using the code provided online
# Change the paths of the input files to rerun the script
rm(list = ls())
library(DistMap)

# Input files #
raw.data = read.csv("Data/Drosophila/dge_raw.txt.gz",sep = "\t",header = F)
rownames(raw.data) = raw.data$V1
raw.data$V1 = NULL
normalized.data = read.csv("Data/Drosophila/dge_normalized.txt.gz", sep = "\t")
insitu.matrix = read.csv("Data/Drosophila/binarized_bdtnp.csv.gz",check.names=F)
geometry = read.csv("Data/Drosophila/geometry.txt.gz",sep = " ")

# Run DistMap
dm = new("DistMap", raw.data=as.matrix(raw.data),
         data = as.matrix(normalized.data),
         insitu.matrix = as.matrix(insitu.matrix),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)

# Save mapped Cells
save(dm, file = "Results_Common/my_mapCells_run.RData")

