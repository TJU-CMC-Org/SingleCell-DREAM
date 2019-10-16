# ************************
# Create folders tree ####
# Default folder locations that are used to save and read files are generated.
# ************************
library(R.utils)
cat("\nCreating folders tree...\n")

# Data folder #
cat("- Data folders\n")
dir.create(path = "Data")
# Nested CV folders
dir.create(path = "Data/CV_folds")

# Results folders LASSO #
cat("- LASSO.topX workflow Results folders...\n")
dir.create(path = "LASSO_topX_workflow/Results_UniquelyMapped_cells_inSituRNAseq")
dir.create(path = "LASSO_topX_workflow/Results_UniquellyMapped_cells_allRNASeq")

# Results for NN will stay in the respective NN folders

# Results folders Common #
cat("- Common Results folders\n\n")
dir.create(path = "Results_Common")
dir.create(path = "Results_Common/Baseline_Method")
dir.create(path = "Results_Common/SubmissionFiles_CV_DistMapTrainTest")
dir.create(path = "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData")

# Download files #
cat("Download data files from: \nhttps://shiny.mdc-berlin.de/DVEX/\n")
download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_raw.txt.gz", 
              destfile = "Data/dge_raw.txt.gz")
download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz",
              destfile = "Data/dge_normalized.txt.gz")
download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/binarized_bdtnp.csv.gz",
              destfile = "Data/binarized_bdtnp.csv.gz")
download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/bdtnp.txt.gz",
              destfile = "Data/bdtnp.txt.gz")
download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/geometry.txt.gz",
              destfile = "Data/geometry.txt.gz")

# Unzip files #
cat("Unzip data files\n")
files2unzip <- list.files(path = "Data/", pattern = ".gz", full.names = TRUE)
sapply(X = files2unzip, FUN = gunzip, remove = FALSE)

cat("\nAll done!\n")
