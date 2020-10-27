#This function is a shortcut to plotting 2D density plots
#using hexbin

plot.hexbin <- function(x, y, xlab, ylab, main){
    require(hexbin)
    require(RColorBrewer)

    if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

    rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

    df <- data.frame(cbind(x, y))
    h <- hexbin(df)

    plot(h, xlab = xlab, ylab = ylab, colramp = rf, main = main)
}