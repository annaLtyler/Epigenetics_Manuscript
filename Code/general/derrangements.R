#generates a matrix of all derrangements for a given number of elements.
#these are all permutations with no fixed points. That is no person gets their
#own hat back.
#perm.mat <- derrangements(8)
#from https://codegolf.stackexchange.com/questions/184943/enumerate-derangements


derrangements <- function(n,y=gtools::permutations(n,n))y[!colSums(t(y)==1:n),]