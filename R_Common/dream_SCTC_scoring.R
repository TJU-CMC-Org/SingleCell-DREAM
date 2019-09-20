#Sets up the enviroment for scoring
#Defines the scoring metrics

library(DistMap)
library(purrr)
library(dplyr)
library(mccr)
library(caret)
library(tidyverse)
# library(synapser)

#Use Attila's code to donwload, load data and initialize the environment for scoring <- syn16782361
initialize <- function(){
  if(!file.exists("init.RData")){
    if(!all(file.exists(c("dge_raw.txt.gz","dge_normalized.txt.gz","binarized_bdtnp.csv.gz","bdtnp.txt.gz","geometry.txt.gz")))){
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_raw.txt.gz",destfile = "dge_raw.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz",destfile = "dge_normalized.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/binarized_bdtnp.csv.gz",destfile = "binarized_bdtnp.csv.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/bdtnp.txt.gz",destfile = "bdtnp.txt.gz")
      download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/geometry.txt.gz",destfile = "geometry.txt.gz")
    }
    
    raw.data = read.table(gzfile("dge_raw.txt.gz",'rt'),
                          sep = "\t",
                          row.names = NULL,
                          stringsAsFactors = F,
                          quote = "")
    raw.data.genes = raw.data$V1
    raw.data$V1 = NULL
    
    raw.data.genes = gsub("'","",raw.data.genes,fixed = T)
    
    raw.data = as.matrix(raw.data)
    rownames(raw.data) = raw.data.genes
    
    normalized.data = read.table(gzfile("dge_normalized.txt.gz",'rt'),
                                 sep = "\t",
                                 row.names = NULL,
                                 stringsAsFactors = F,
                                 quote = "")
    
    normalized.data.genes = normalized.data$row.names
    normalized.data$row.names = NULL
    
    normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)
    
    normalized.data = as.matrix(normalized.data)
    rownames(normalized.data) = normalized.data.genes
    
    stopifnot(all(normalized.data.genes == raw.data.genes))
    
    insitu.matrix = read.table(gzfile("binarized_bdtnp.csv.gz",'rt'), sep = ",",header = T)
    
    insitu.genes_orig <- colnames(insitu.matrix)
    
    #this is not needed for the normalized data
    insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
    insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)
    
    insitu.matrix = as.matrix(insitu.matrix)
    colnames(insitu.matrix) = insitu.genes
    
    stopifnot(all(insitu.genes %in% raw.data.genes))
    
    geometry = read.csv(gzfile("geometry.txt.gz",'rt'),sep = " ")
    
    colnames(geometry) = c("x","y","z")
    
    #close gz properly
    #closeAllConnections()
    
    dm <<- new("DistMap",
             raw.data=raw.data,
             data=normalized.data,
             insitu.matrix=insitu.matrix,
             geometry=as.matrix(geometry))
    
    dm <<- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
    
    dm <<- mapCells(dm)
    #Thank you Attila!
    
    #GROUND TRUTH
    
    ground.truth <<- t(apply(dm@mcc.scores,2,order,decreasing=TRUE))[,1:10]
    ambig.locations <<- t(apply(dm@mcc.scores,2,sort,decreasing=TRUE))[,1:2]
    ambig.locations <<- which(ambig.locations[,1] == ambig.locations[,2])
    
    
    #map every cell to its d84 value
    d84 <<- seq(nrow(ground.truth)) %>% map_dbl(function(j){
      #map every position to the norm of the difference in the geometry and calculate the mean
      ground.truth[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
    })
    
    save(dm,ground.truth,d84,ambig.locations,file="init.RData")
    
  } else {
    load("init.RData",envir = .GlobalEnv)
  }
}

#Scoring function. 
#Input: path to the results .csv file (character), number of subchallenge (integer)
#Output: vector of scores (s1,s2,s3)
score <- function(path,sub){
  
  if (!exists("dm")) initialize()
  
  submission <- read.csv(path,header=FALSE,stringsAsFactors = FALSE)
  
  #separate the gene names from the location predictions
  gene.lines <- (4-sub)*2
  genes <- submission %>% slice(1:gene.lines)
  locations <- submission %>% slice(-1:-gene.lines)
  
  #preprocess genes and locations, remove NAs, sort locations by cellid
  genes <- genes %>% select(-1) %>% unlist %>% as.character
  locations <- locations[order(locations[,1]),] %>% select(-1) %>% apply(2,as.numeric)
  
  genes = gsub("-",".",genes,fixed = T)
  genes = gsub("(spl)",".spl.",genes,fixed = T)
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations)) %>% map_dbl(function(j){
    vals <- locations[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84/dsub
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
  #get binarized data from distmap
  ts <- data.frame(t(dm@binarized.data))
  #select binarized data only for the submitted subset of genes
  reduced.ts <- ts %>% select(genes)
  
  #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],reduced.insitu[locations[.x,1],]))
  
  #do not take into account the cells with ambiguous locations
  s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.locations])
  
  #s2
  #do not take into account the cells with ambiguous locations
  s2<- mean(pk[-ambig.locations])
  
  #s3
  
  #comparing rnaseq and fluorescence data using true locations
  true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[-ambig.locations,1],.x],reduced.ts[-ambig.locations,.x]))
  #.. using submitted locations
  competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.locations,1],.x],reduced.ts[-ambig.locations,.x]))
  
  #do not take into account the cells with ambiguous locations
  s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
  
  return(c(s1,s2,s3))
}

#Scoring function with bootstraping
#Input: path to the results .csv file (character), number of subchallenge (integer), number of bootstraps (integer, optional)
#Output: data frame with scores
score.bootstrapped <- function(path,sub,nboot=1000){
  
  if (!exists("dm")) initialize()
  
  submission <- read.csv(path,header=FALSE,stringsAsFactors = FALSE)
  
  #separate the gene names from the location predictions
  gene.lines <- (4-sub)*2
  genes <- submission %>% slice(1:gene.lines)
  locations <- submission %>% slice(-1:-gene.lines)
  
  #preprocess genes and locations, remove NAs, sort locations by cellid
  genes <- genes %>% select(-1) %>% unlist %>% as.character
  locations <- locations[order(locations[,1]),] %>% select(-1) %>% apply(2,as.numeric)
  
  #fix incompatibility
  genes = gsub("-",".",genes,fixed = T)
  genes = gsub("(spl)",".spl.",genes,fixed = T)
  
  #remove ambiguous locations
  locations.n <- locations[-ambig.locations,]
  ground.truth.n <- ground.truth[-ambig.locations,]
  
  #do the same mapping for the submission as for d48
  dsub <- seq(nrow(locations.n)) %>% map_dbl(function(j){
    vals <- locations.n[j,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth.n[j,1],])^2))) %>% mean
  })
  
  #calculate relative precision
  pk <- d84[-ambig.locations]/dsub
  
  #s1
  
  #select fluorescence data only for the submitted subset of genes
  reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
  #get binarized data from distmap and remove abiguous locations
  ts <- data.frame(t(dm@binarized.data))[-ambig.locations,]
  #select binarized data only for the submitted subset of genes
  reduced.ts <- ts %>% select(genes)
  
  #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
  mccrs <- seq(nrow(locations.n)) %>% map_dbl(~mccr(reduced.insitu[ground.truth.n[.x,1],],reduced.insitu[locations.n[.x,1],]))
  
  #bootstrapping
  samples <- seq(nboot) %>% map_dfr(function(seed){
    set.seed(seed)
    bootstrap <- sample.int(nrow(locations.n),replace=TRUE)
    #s1
    s1.b <- sum((pk[bootstrap]/sum(pk[bootstrap])) * mccrs[bootstrap])
    #s2
    s2.b <- mean(pk[bootstrap])
    
    #s3
    #since we bootstrap by locations, these must be recalculated
    true.mccs.b <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth.n[bootstrap,1],.x],reduced.ts[bootstrap,.x]))
    #submitted locations
    competitor.mccs.b <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations.n[bootstrap,1],.x],reduced.ts[bootstrap,.x]))
    
    #here i assumed that the denominator is the sum of true mccs
    s3.b <- sum((true.mccs.b/sum(true.mccs.b)) * competitor.mccs.b)
    
    
    data.frame(s1=s1.b,s2=s2.b,s3=s3.b)
  })
  
  return(samples)
}

#wrapper for summarising the bootstrapped scores
score.bootstrapped.summary <- function(path,sub,nboot=1000){
  score.bootstrapped(path,sub,nboot) %>% summarise(mean(s1), sd(s1), mean(s2), sd(s2), mean(s3), sd(s3)) %>% as.numeric
}

#computes the Bayes factor between two results using bootstrapped scores
bayes.bootstrap <- function(path1,path2,sub,nboot=1000){
  
  samples1 <- score.bootstrapped(path1,sub,nboot)
  samples2 <- score.bootstrapped(path2,sub,nboot)
  
  wins <- colSums(samples1 >= samples2)
  B <- wins/(nboot-wins)
  
  return(B)
}

#input: a csv file with a sid column and a team column, requires synapse login
bootstraped.ranks <- function(submissions,sub){
  
  s <- read.csv(submissions,stringsAsFactors = F)
  s$sid <- as.character(s$sid)
  
  #if not logged in synapse, load in the variable files the paths to the submissions
  files <- s$sid %>% map_chr(~synGetSubmission(.x)$filePath)
  
  #evaluate
  eval.boot <-files %>% map(~score.bootstrapped(.x,sub))
  
  #rank on each score separately and reduce to sum
  #need to make scores negative in order to properly rank them
  ranks <- seq(3) %>% map(function(score){
    ranks<- eval.boot %>% map_dfc(~-.x[,score]) %>% apply(1,rank) %>% t
    colnames(ranks) <- s$team  
    return(ranks)
  }) %>% reduce(`+`)
  
  ranks <- (ranks/3) %>% apply(1,rank) %>% t
  save(ranks,file=paste0("sc",sub,"ranks.Rdata"))
  
  #draw the boxplot
  avg.ranks <- ranks %>% colMeans %>% rank  
  ordering <- order(avg.ranks)
  pdf(paste0("sc",sub,"_final_boxplot.pdf"), width=11, height=8)
  par(mar=c(5,10.5,4,2) + 0.1)
  boxplot(ranks[,ordering],horizontal=T,las=2,at=rev(1:ncol(ranks)), xlab="Rank")
  
  factors <- map2_dfr(ordering[-length(ordering)],ordering[-1],function(c1,c2){
    win <- sum(ranks[,c1] < ranks[,c2])
    lose <- sum(ranks[,c2] < ranks[,c1])
    BF <- win/lose
    data.frame(c1=c1,c2=c2,BF=BF)
  })
  
  abh <- ncol(ranks) - factors$c1[which(factors$BF>=3)[1]] + 0.5
  abline(h=abh,lwd=2)
  
  dev.off()
  
  #report final ranking
  result <- mutate(s, rank=avg.ranks) %>% arrange(rank)
  write.csv(result,file=paste0("sc",sub,"_final_table.csv"),row.names=F)
  
  return(mutate(s, rank=avg.ranks) %>% arrange(rank))
}

#score function for the 10 fold crossvalidation
#the pattern should be provided in glue format, where the fold number is represented by {.} or {.x}
score.post.folds <- function(pattern, sub){
    if (!exists("dm")) initialize()
    
    #read all folds
    all.submissions <- seq(10) %>% map(~read.csv(str_glue(pattern),header=FALSE,stringsAsFactors = FALSE))
    
    #for each fold get scores
    all.submissions %>% map_dfr(function(submission){
        #separate the gene names from the location predictions
        gene.lines <- (4-sub)*2
        genes <- submission %>% slice(1:gene.lines)
        locations <- submission %>% slice(-1:-gene.lines)
        #BCBU hack
        #locations <- locations[which(locations[,2] != ""),]
        
        #preprocess genes and locations, remove NAs, sort locations by cellid
        #mutations cause cancer, but anyway ...
        genes <- genes %>% select(-1) %>% unlist %>% as.character
        
        #fix incompatibility
        genes = gsub("-",".",genes,fixed = T)
        genes = gsub("(spl)",".spl.",genes,fixed = T)
        
        locations <- locations[order(locations[,1]),] 
        cell.ids <- locations[,1]
        locations <- locations %>% select(-1) %>% apply(2,as.numeric)
        
        #select from d84 the falues from the current test fold 
        d84.fold <- d84[cell.ids]
        
        #do the same mapping for the submission as for d84
        dsub <- cell.ids %>% imap_dbl(function(j,iter){
            vals <- locations[iter,] %>% map_dbl(~sqrt(sum((dm@geometry[.x,] - dm@geometry[ground.truth[j,1],])^2))) %>% mean
        })
        
        #calculate relative precision
        pk <- d84.fold/dsub
        
        #s1
        
        #select fluorescence data only for the submitted subset of genes
        reduced.insitu <- data.frame(dm@insitu.matrix) %>% select(genes)
        #get binarized data from distmap
        ts <- data.frame(t(dm@binarized.data))
        #select binarized data only for the submitted subset of genes
        reduced.ts <- ts %>% select(genes)
        
        
        #map every cell location prediction to the MCC between the ground truth location and the predicted most likely position, using the submitted subset of genes 
        mccrs <- cell.ids %>% imap_dbl(~mccr(reduced.insitu[ground.truth[.x,1],],reduced.insitu[locations[.y,1],]))
        
        ambig.indexes <- which(cell.ids %in% ambig.locations)
        
        #do not take into account the cells with ambiguous locations
        s1 <- sum(((pk/sum(pk)) * mccrs)[-ambig.indexes])
        
        #s2
        #do not take into account the cells with ambiguous locations
        s2 <- mean(pk[-ambig.indexes])
        
        #s3
        
        #comparing rnaseq and fluorescence data using true locations
        true.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[ground.truth[cell.ids[-ambig.indexes],1],.x],reduced.ts[cell.ids[-ambig.indexes],.x]))
        #.. using submitted locations
        competitor.mccs <- seq(ncol(reduced.ts)) %>% map_dbl(~mccr(reduced.insitu[locations[-ambig.indexes,1],.x],reduced.ts[cell.ids[-ambig.indexes],.x]))
        
        #do not take into account the cells with ambiguous locations
        s3 <- sum(((true.mccs/sum(true.mccs)) * competitor.mccs))
        
        data.frame(s1=s1,s2=s2,s3=s3)
        
    })
    
}



