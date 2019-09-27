rm(list = ls())
library(DistMap)
raw.data = read.csv("../../Data/dge_raw.txt",sep = "\t",header = F)
rownames(raw.data) = raw.data$V1
raw.data$V1 = NULL

normalized.data = read.csv("../../Data/dge_normalized.txt", sep = "\t")
insitu.matrix = read.csv("../../Data/binarized_bdtnp.csv",check.names=F)

geometry = read.csv("../../Data/geometry.txt",sep = " ")
dm = new("DistMap", raw.data=as.matrix(raw.data),
         data = as.matrix(normalized.data),
         insitu.matrix = as.matrix(insitu.matrix),
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))

dm <- mapCells(dm)


# Extract locations for RNAseq data
mccScores <- dm@mcc.scores
colnames(mccScores) <- colnames(normalized.data)
dim(mccScores)

# Index of cells position RNAseq position
index.RNAseq <- apply(X = mccScores, MARGIN = 2, FUN = function(coli){
  thresh_i = max(coli) * 0.95
  return(which(coli >= thresh_i))
})

trainmatrix = matrix(, nrow = length (unlist (index.RNAseq)), ncol = (length (rownames(normalized.data)) + dim (geometry)[2]))
colnames (trainmatrix) = c (colnames (geometry), rownames(normalized.data));
addrow = 1;
for (cell_i in 1:length (index.RNAseq)) # iterate over each cell
{
  for (position_i in 1:length(index.RNAseq[[cell_i]]))  # iterate over each position
  {
      rowdata_geometry = geometry[index.RNAseq[[cell_i]][position_i],];
      rowdata_genedata = normalized.data[,cell_i];
      
      trainmatrix[addrow,] = as.numeric (c(rowdata_geometry, rowdata_genedata));
      addrow = addrow + 1
  }
}

write.table(trainmatrix, file="trainingTable.txt", sep="\t", row.names=FALSE, col.names=TRUE)
