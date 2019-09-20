#################################################
# SUMMARY
# This script runs baseline method
# Baseline method randomly selects 20 / 40 / 60 genes.
#################################################
rm(list = ls())

# ***************
# USER INPUT ####
# ***************
numGenes <- c(20, 40, 60)
numFolds <- 10

path2save <- "Results_Common/Baseline_Method/"

# *****************
# CALCULATIONS ####
# *****************
library(data.table)


# Import raw and normalised data decompress the downloaded dataset and import to data.frame.
# Important: quote is set to empty string, this way we can prevent the mis-interpretation of a gene name, that contains an apostrophe

# We convert the data.frame to matrix: each column represents a cell and each row corresponds to a gene.
raw.data = read.table("Data/dge_raw.txt", sep = "\t", row.names = NULL,
                      stringsAsFactors = F, quote = "")
raw.data.genes = raw.data$V1
raw.data$V1 = NULL

# Let’s fix the gene names that contains apostrophe – this would generate issues
# gene names with apostrophes
# print(grep("'",raw.data.genes,value = T,fixed = T))
# [1] "beta'COP" "PP2A-B'"
raw.data.genes = gsub("'","",raw.data.genes,fixed = T)

raw.data = as.matrix(raw.data)
rownames(raw.data) = raw.data.genes

# Repeat for the normalised data
normalized.data = read.table("Data/dge_normalized.txt", sep = "\t",
                             row.names = NULL, stringsAsFactors = F, quote = "")
normalized.data.genes = normalized.data$row.names
normalized.data$row.names = NULL

# gene names with apostrophes
# print(grep("'",normalized.data.genes,value = T,fixed = T))

normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)

normalized.data = as.matrix(normalized.data)
rownames(normalized.data) = normalized.data.genes

# MY ASSUMPTION #
colnames(raw.data) <- colnames(normalized.data)

# Check that the gene names are identical in the raw and normalised dataset
stopifnot(all(normalized.data.genes == raw.data.genes))

# Import in situ datasets
insitu.matrix = read.table("Data/binarized_bdtnp.csv", sep = ",",header = T)

## Warning in read.table(gzfile("binarized_bdtnp.csv.gz", "rt"), sep = ",", :
## seek on a gzfile connection returned an internal error

## Warning in read.table(gzfile("binarized_bdtnp.csv.gz", "rt"), sep = ",", :
## seek on a gzfile connection returned an internal error

insitu.genes_orig <- colnames(insitu.matrix)

# Match the gene names across datasets

# The following few lines of code checks if there is any mismatch between the gene names.

# We will find that 2 gene names were changed during the import of the dataset and the - (dash) character was changed to .(dot).
# 
# Further, one gene names contained brackets (), which were replaced by dashes .

# 2 gene names are not matched:
missingGenes = insitu.genes_orig[which(!insitu.genes_orig %in% normalized.data.genes)]
# print(missingGenes)

## [1] "Blimp.1"      "E.spl.m5.HLH"

# this was reported by Nikos
# lets fix this by changing the . characters in the gene names to -
insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
# also replace .spl. --> (spl)
insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)

# assert that all institu genes appear in the gene names
stopifnot(all(insitu.genes %in% raw.data.genes))

# Now we can rename the genes in the institu.matrix with the correct names:

insitu.matrix = as.matrix(insitu.matrix)
colnames(insitu.matrix) = insitu.genes

# Select features randomly
featuresInfo <- data.table(numGenes = rep(numGenes, each = numFolds), 
                           fold = 1:numFolds)

for(i in 1:nrow(featuresInfo)){
    importantGenes <- sample(x = colnames(insitu.matrix), 
                             size = featuresInfo$numGenes[i])
    write.table(x = importantGenes, 
                file = paste0(path2save, 
                              featuresInfo$numGenes[i], "genes", 
                              "_CV_", featuresInfo$fold[i], ".txt"), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
}
















