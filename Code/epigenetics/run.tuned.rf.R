#This function takes in the features for a Random forest
#and the output of tune.maxnodes and runs the random forest
#at the best number of maxnodes.


run.tuned.rf <- function(y, x, tuned.param, ntree = 500){

    min.error <- which.min(tuned.param[,1])
    max.exp <- which.max(tuned.param[,2])
    max.dist <- which.max(tuned.param[,2] - tuned.param[,1])
    max.nodes <- round(mean(c(min.error, max.exp, max.dist)))

    df <- data.frame(x)

    rf.result <- randomForest(y~., data = df, maxnodes = max.nodes, ntree = ntree, localImp = TRUE)
    return(rf.result)
}