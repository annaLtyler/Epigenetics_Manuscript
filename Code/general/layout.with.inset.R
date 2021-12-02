#This function generates a layout matrix for a 
#plot with an inset, where the inset's coordinates
#are determined by x and y coordinates.

layout.with.inset <- function(nrow = 10, ncol = 10, inset.x.min = 3, inset.x.max = 5, 
inset.y.min = 3, inset.y.max = 5){

    layout.mat <- matrix(rep(1, (nrow*ncol)), nrow = nrow, ncol = ncol)

    #find the matrix entries for the inset
    x.coord <- inset.x.min : inset.x.max
    #y.coord <- (nrow - inset.y.max + 1) : (nrow - inset.y.min + 1)
    y.coord <- ((nrow - inset.y.max) - 1): ((nrow - inset.y.min) + 1)
    inset.coord <- cbind(rep(x.coord, length = length(y.coord)), rep(y.coord, each = length(x.coord)))
    for(i in 1:nrow(inset.coord)){
        layout.mat[inset.coord[i,2], inset.coord[i,1]] <- 2
    }

    #layout(layout.mat)
    #layout.show(2)

    return(layout.mat)
}