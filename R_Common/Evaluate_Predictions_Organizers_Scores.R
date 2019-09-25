rm(list = ls())
library(foreach)
library(doMC)
library(ggfortify)
library(data.table)
source("R_Common/dream_SCTC_scoring.R")


# ***************
# USER INPUT ####
# ***************
# Folder to save results
folder2save <- "Results_Common/"


# *****************
# CALCULATIONS ####
# *****************

# Generate scores CV ####
scoresCalculator <- function(subch_i, pattern_files){
    scores_sub_i <- score.post.folds(pattern = pattern_files[subch_i], sub = subch_i)
    scores_sub_i <- data.table(scores_sub_i)
    scores_sub_i$Subchallenge <- subch_i
    scores_sub_i$Method <- paste(strsplit(x = basename(pattern_files[subch_i]), 
                                          split = ".", fixed = TRUE)[[1]][1], 
                                 strsplit(x = strsplit(x = pattern_files[subch_i], 
                                          split = "/", fixed = TRUE)[[1]][2], 
                                          split = "_", fixed = TRUE)[[1]][3], sep = ".")
    return(scores_sub_i)
}

# @ TrainTest ####
patterns_TrainTest_LASSO <- c("Results_Common/SubmissionFiles_CV_DistMapTrainTest/Lasso.DistMap.60.CV.{.}.csv", 
                              "Results_Common/SubmissionFiles_CV_DistMapTrainTest/Lasso.DistMap.40.CV.{.}.csv",
                              "Results_Common/SubmissionFiles_CV_DistMapTrainTest/Lasso.DistMap.20.CV.{.}.csv")

system.time({
    allScoresTrainTest_LASSO <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                               pattern_files = patterns_TrainTest_LASSO))
})

patterns_TrainTest_NN <- c("Results_Common/SubmissionFiles_CV_DistMapTrainTest/NeuralNets.DistMap.60.CV.{.}.csv",
                           "Results_Common/SubmissionFiles_CV_DistMapTrainTest/NeuralNets.DistMap.40.CV.{.}.csv",
                           "Results_Common/SubmissionFiles_CV_DistMapTrainTest/NeuralNets.DistMap.20.CV.{.}.csv")


system.time({
    allScoresTrainTest_NN <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                                     pattern_files = patterns_TrainTest_NN))
})


patterns_TrainTest_Random <- c("Results_Common/SubmissionFiles_CV_DistMapTrainTest/Random.DistMap.60.CV.{.}.csv", 
                               "Results_Common/SubmissionFiles_CV_DistMapTrainTest/Random.DistMap.40.CV.{.}.csv",
                               "Results_Common/SubmissionFiles_CV_DistMapTrainTest/Random.DistMap.20.CV.{.}.csv")

system.time({
    allScoresTrainTest_Random <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                                  pattern_files = patterns_TrainTest_Random))
})



# @ Provided Binary table ####
patterns_UsingProvidedBinaryData_LASSO <- 
    c("Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Lasso.DistMap.60.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Lasso.DistMap.40.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Lasso.DistMap.20.CV.{.}.csv")

system.time({
    allScoresProvidedBinaryData_LASSO <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                                     pattern_files = patterns_UsingProvidedBinaryData_LASSO))
})


patterns_UsingProvidedBinaryData_NN <- 
    c("Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/NeuralNets.DistMap.60.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/NeuralNets.DistMap.40.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/NeuralNets.DistMap.20.CV.{.}.csv")

system.time({
    allScoresProvidedBinaryData_NeuralNets <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                                              pattern_files = patterns_UsingProvidedBinaryData_NN))
})


patterns_UsingProvidedBinaryData_Random <- 
    c("Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Random.DistMap.60.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Random.DistMap.40.CV.{.}.csv", 
      "Results_Common/SubmissionFiles_CV_DistMapOnTestCells_UsingProvidedBinaryData/Random.DistMap.20.CV.{.}.csv")

system.time({
    allScoresProvidedBinaryData_Random <- rbindlist(l = lapply(X = c(1:3), FUN = scoresCalculator, 
                                                               pattern_files = patterns_UsingProvidedBinaryData_Random))
})


# Combine all ####
allScores <- rbindlist(l = list(allScoresTrainTest_LASSO, 
                                allScoresTrainTest_NN, 
                                allScoresTrainTest_Random, 
                                allScoresProvidedBinaryData_LASSO, 
                                allScoresProvidedBinaryData_NeuralNets, 
                                allScoresProvidedBinaryData_Random))

allScores[, Method := gsub(pattern = ".DistMapTrainTest", replacement = "_testCell_selGenes", x = Method)]
allScores[, Method := gsub(pattern = ".DistMapOnTestCells", replacement = "_testCell_PB", x = Method)]


# Comparison plots ####
allScores <- melt.data.table(data = allScores, 
                             id.vars = c("Subchallenge", "Method"), 
                             variable.name = "Metric", value.name = "Score")

allScores[, meanValue := mean(Score), by = c("Subchallenge", "Method", "Metric")]

myplots <- foreach(sub_i = unique(allScores$Subchallenge)) %do% {
    ggplot(data = allScores[Subchallenge == sub_i],
           mapping = aes(x = Metric, y = Score, fill = Method)) + 
        geom_boxplot() +
        ylab(label = "Score (Nested CV)") +
        # scale_y_continuous(breaks = seq(from = 0.5, to = 1.2, by = 0.1), limits = c(0.5, 1.2)) +
        # geom_point() +
        ggtitle(paste("Subchallenge =", sub_i)) + 
        theme_bw()
}
print(autoplot(object = myplots))

ggsave(filename = paste0(folder2save, "MethodsComparison_NestedCV_Distribution_OrganizersScores.jpg"), 
       plot = autoplot(object = myplots), width = 8, height = 10)


myplotsMean <- foreach(sub_i = unique(allScores$Subchallenge)) %do% {
    ggplot(data = allScores[Subchallenge == sub_i],
           mapping = aes(x = Metric, y = meanValue, fill = Method)) + 
        geom_bar(stat = "identity", position = position_dodge()) +
        ylab(label = "Mean score (Nested CV)") + 
        scale_y_continuous(breaks = seq(from = 0, to = 1.1, by = 0.2), limits = c(0, 1.1)) +
        # geom_point() +
        ggtitle(paste("Subchallenge =", sub_i)) +
        theme_bw() #+ theme(legend.position = "bottom")
}

print(autoplot(object = myplotsMean))
ggsave(filename = paste0(folder2save, "MethodsComparison_NestedCV_Mean_OrganizersScores.jpg"), 
       plot = autoplot(object = myplotsMean), width = 7, height = 10)

fwrite(x = allScores, file = paste0(folder2save, "OrganizersScoreFunctionsResults.csv"), 
       quote = FALSE)
