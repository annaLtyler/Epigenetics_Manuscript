#This function runs a random forest by sampling genes
#the given number of times. gene.id is a vector that 
#corresponds to the rows of the feature matrix. Rows
#from the same gene are labeled with the same index
#so that all rows from a single gene are sampled together.

run.rf <- function(y, feature.mat, trials = 10, n.genes = 100, gene.id, 
ntree = 1000, maxnodes = ncol(feature.mat)*3, file.label = "Random.Forest"){

total.genes <- length(unique(gene.id))
all.rf <- vector(mode = "list", length = trials)
for(tr in 1:trials){
  forest.file <- file.path("results", paste0(file.label, ".", tr, ".RDS"))
  
  if(!file.exists(forest.file)){
    sampled.genes <- sample(1:total.genes, n.genes)
    sub.idx <- which(gene.id %in% sampled.genes)
    sub.feature.mat <- feature.mat[sub.idx,]
    sub.y <- y[sub.idx,]

    has.na <- as.logical(length(which(is.na(sub.feature.mat))))

    if(has.na){
      imp.file <- file.path("results", paste0(file.label, ".Imputed.", tr, ".RDS"))
      if(!file.exists(imp.file)){
          imp <- rfImpute(sub.feature.mat, sub.y)
          saveRDS(imp, imp.file)
      }else{
          imp <- readRDS(imp.file)
      }
    }else{
      imp <- cbind(sub.y, sub.feature.mat)
    }

    #tune the mtry variable
    #test <- tuneRF(imp[,2:ncol(imp)], imp[,1], stepFactor = 2, trace = FALSE, 
    #plot = FALSE, doBest = FALSE)
    #mtry <- test[which.min(test[,2]),1]

    imp.df <- data.frame(imp)
    
    if(is.null(ntree) && is.null(maxnodes)){
      #rf <- randomForest(imp[,2:ncol(imp)], imp[,1], mtry = mtry, localImp = TRUE)  
      #rf <- randomForest(imp[,1]~imp[,2:ncol(imp)], mtry = mtry, localImp = TRUE)  

      #using the formula format allows us to look at interactions
      #using randomForestExplainer
      rf <- randomForest(sub.y~., data = imp.df, localImp = TRUE)  
    }else{
      #rf <- randomForest(imp[,1]~imp[,2:ncol(imp)], mtry = mtry, localImp = TRUE,
      #ntree = ntree, maxnodes = maxnodes)
      #rf <- randomForest(imp[,2:ncol(imp)], imp[,1], mtry = mtry, localImp = TRUE,
      #ntree = ntree, maxnodes = maxnodes)

      #using the formula format allows us to look at interactions
      #using randomForestExplainer
      rf <- randomForest(sub.y~., data = imp.df, localImp = TRUE, ntree = ntree,
      maxnodes = maxnodes)  
    }
    
    #varImpPlot(rf)
    #compare.importance(rf)
    saveRDS(rf, forest.file)
    all.rf[[tr]] <- rf
  }else{
    rf <- readRDS(forest.file)
    all.rf[[tr]] <- rf
  }
}

if(trials == 1){
  return(all.rf[[1]])
}else{
  return(all.rf)
}

}