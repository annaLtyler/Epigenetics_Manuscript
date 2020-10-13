#this function plots correlations between chromatin state proportion
#and inbred expression 


state.prop.v.expression <- function(chrom.state.prop, centered.expr, name.color.key, start.feature, end.feature, upstream.buffer, downstream.buffer){
	
	num.strains <- ncol(chrom.state.prop[[1]])
	num.states <- nrow(chrom.state.prop[[1]])
	# if(relative.chrom.state){
		# relative.text <- "relative.to.mean"
		# }else{
		# relative.text <- ""	
		# }
		
	get.strain.stats <- function(ind.chrom.states, strain.chrom.locale){
		if(!is.na(ind.chrom.states[1])){
			return(ind.chrom.states[,strain.chrom.locale])
			}else{
			return(rep(NA, num.states))
			}
		}
	
	pdf(paste("State.Expression.Comparisons.", start.feature, ".", end.feature, ".", upstream.buffer, ".", downstream.buffer, ".pdf", sep = ""), width = 9, height = 6)

	all.models <- vector(mode = "list", length = nrow(name.color.key))
	names(all.models) <- name.color.key[,4]
	for(ms in 1:nrow(name.color.key)){
		strain.models <- vector(mode = "list", length = num.states)
		names(strain.models) <- state.key[,2]
		if(!is.na(name.color.key[ms,5])){
			strain.chrom.locale <- which(colnames(chrom.state.prop[[1]]) == name.color.key[ms,4])
			strain.state.mat <- matrix(unlist(lapply(chrom.state.prop, function(x) get.strain.stats(x, strain.chrom.locale))), ncol = num.states, byrow = TRUE)
			colnames(strain.state.mat) <- 1:num.states
			rownames(strain.state.mat) <- names(chrom.state.prop)
			
			str.expr.locale <- which(names(centered.expr[[1]]) == name.color.key[ms,4])
			str.expr <- unlist(lapply(centered.expr, function(x) x[str.expr.locale]))
			
			layout.mat <- get.layout.mat(num.states)
			layout(layout.mat)
			for(st in 1:num.states){
				plot(strain.state.mat[,st], str.expr, xlab = "Proportion State", ylab = "Expression", main = paste(state.key[st,2],"\n", state.key[st,3], sep = ""))
				model <- lm(str.expr~strain.state.mat[,st])
				abline(model)
				strain.models[[st]] <- model
				}
			all.models[[ms]] <- strain.models
			mtext(paste("Chromatin States v. Expression:", name.color.key[ms,1]), outer = TRUE, line = -2)
			
			# strain.coef.locale <- which(names(all.cis.coef[[1]]) == name.color.key[ms,5])
			# strain.coef <- unlist(lapply(all.cis.coef, function(x) x[[strain.coef.locale]]))
			# common.genes <- intersect(rownames(strain.state.mat), names(strain.coef))
			# gene.order1 <- match(common.genes, rownames(strain.state.mat))
			# gene.order2 <- match(common.genes, names(strain.coef))
			# gene.table <- cbind(rownames(strain.state.mat)[gene.order1], names(strain.coef)[gene.order2])
			# # identical(gene.table[,1], gene.table[,2])	
			
			# layout(layout.mat)
			# for(st in 1:num.states){
				# plot(strain.state.mat[gene.order1,st], strain.coef[gene.order2], xlab = "Proportion State", ylab = "DO Coefficient", main = state.key[st,2])
				# # smoothScatter(strain.state.mat[gene.order1,st], strain.coef[gene.order2], xlab = "Proportion State", ylab = "DO Coefficient", main = state.key[st,2])
				# }
			# mtext(paste("Chromatin States v. DO Coefficients:", name.color.key[ms,2]), outer = TRUE, line = -1.5)
			}
		}
	dev.off()
	
	pdf(paste("State.Expression.Comparisons.Models.", start.feature, ".", end.feature, ".", upstream.buffer, ".", downstream.buffer, ".pdf", sep = ""), width = 5, height = 6)
	layout(matrix(c(1,2), nrow = 2), heights = c(5/6, 1-5/6))
	for(i in 1:length(all.models[[1]])){
		par(mar = c(4,4,4,4))
		plot.new()
		plot.window(xlim = c(-10, 10), ylim = c(-10, 10))
		for(j in 1:length(all.models)){
			abline(all.models[[j]][[i]], col = name.color.key[j,3])
			}
		mtext(names(all.models[[1]])[i], cex = 2)
		mtext("State Amount (Arbitraty Units)", side = 1, line = 2.5)
		mtext("Expression (Arbitraty Units)", side = 2, line = 2.5)
		axis(1);axis(2)
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		legend(x = 0, y = 1, legend = name.color.key[,1], col = name.color.key[,3], lty = 1, horiz = TRUE, cex = 0.5, lwd = 2)
		}
	dev.off()


}