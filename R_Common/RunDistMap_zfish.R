# Rerunning distMap using the code provided online
# Change the paths of the input files to rerun the script
rm(list = ls())
library(DistMap)

# Input files #
raw.data = read.table(file = "Data/zfish/data/zfish.raw.tsv", sep = "\t", header = TRUE)
normalized.data = read.table("Data/zfish/data/zfish.normalized.tsv", sep = "\t")
insitu.matrix = read.table("Data/zfish/data/insitus.tsv", sep = "\t", header = TRUE)
geometry = read.csv("Data/zfish/data/pseudo.geometry.tsv",sep = "\t", header = TRUE)

# Run DistMap
dm = new("DistMap", raw.data=as.matrix(raw.data),
         data = as.matrix(normalized.data),
         insitu.matrix = as.matrix(insitu.matrix),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)

# Save mapped Cells
save(dm, file = "Results_Common/my_mapCells_run.RData")


library(data.table)

