compress.window.cor <- function(window.cor){

    window.x <- unlist(lapply(window.cor, function(x) as.numeric(colnames(x))))
    window.beta <- unlist(lapply(window.cor, function(x) x[1,]))
    #window.sd <- unlist(lapply(window.cor, function(x) x[2,]))
    #window.t <- unlist(lapply(window.cor, function(x) x[3,]))
        
    x.order <- order(window.x)

    full.curve <- cbind(window.x[x.order], window.beta[x.order])
    #plot(full.curve[,1], full.curve[,2], type = "l")
    df <- data.frame(full.curve)
    colnames(df) <- c("coord", "val")
    condensed.curve <- aggregate(val~coord, df, mean)
    #plot(condensed.curve[,1], condensed.curve[,2], type = "l")

    vals <- condensed.curve[,2]
    names(vals) <- condensed.curve[,1]

    #plot(u_pts, u_curve, type = "l")
    return(vals)
}