#calculate correlation between DO coefficients
#and inbred mean expression

expr.coef.cor <- function(do.eqtl, expression.list, coef.list, lod.list = c(10, 20, 50), name.color.key){

	max.lod <- unlist(lapply(do.eqtl, function(x) max(x[[3]])))
	names(max.lod) <- names(do.eqtl)
	
	for(l in 1:length(lod.list)){
		
		good.lod.locale <- which(max.lod >= lod.list[l])
		all.cor <- rep(NA, length(good.lod.locale))
		cat("\n")
		for(i in 1:length(good.lod.locale)){
			report.progress(i, length(good.lod.locale))
			gene.id <- do.eqtl[[good.lod.locale[i]]][[1]][1,1]
			
			gene.expr.locale <- which(names(expression.list) == gene.id)
			
			if(length(gene.expr.locale) > 0){
				gene.expr <- expression.list[[gene.expr.locale]]
				
				if(!is.na(gene.expr[1])){
					gene.coef.locale <- which(names(coef.list) == gene.id)
					
					if(length(gene.coef.locale) > 0){
						gene.coef <- coef.list[[gene.coef.locale]]	
						
						if(var(gene.coef) > 0 && var(gene.expr) > 0){
							expr.order <- match(name.color.key[,4], names(gene.expr))
							coef.order <- match(name.color.key[,5], names(gene.coef))
													
							# plot(gene.coef[coef.order], gene.expr[expr.order], col = col.table[,3], cex = 2, pch = 16)
							model <- lm(gene.expr[expr.order]~gene.coef[coef.order])
							# abline(model)
							all.cor[i] <- summary(model)$adj.r.squared
							}
						}
					}
				}
			}
	
		cat("plotting to", paste("expression_v_DO_coef_by_gene.LOD.", lod.list[l], ".pdf", sep = ""), "\n")
		pdf(paste("expression_v_DO_coef_by_gene.LOD.", lod.list[l], ".pdf", sep = ""), width = 4, height = 4)
		hist(all.cor, breaks = 100, main = "", xlab = "Correlation", xlim = c(-0.2, 1))
		dev.off()
			
		centered.expr <- lapply(expression.list, function(x) x - mean(x))
			
		common.genes <- intersect(names(centered.expr)[good.lod.locale], names(coef.list))
		gene.order1 <- match(common.genes, names(centered.expr))
		gene.order2 <- match(common.genes, names(coef.list))
		# identical(names(expression.list)[gene.order1], names(coef.list)[gene.order2])
		layout.mat <- get.layout.mat(length(group.mean.expr[[1]]))
		jpeg(filename = paste("strain_expression_v_DO_coef.LOD.", lod.list[l],".jpg", sep = ""), width = ncol(layout.mat)*250, height = nrow(layout.mat)*250)
		layout(layout.mat)
		for(i in 1:nrow(name.color.key)){
			if(!is.na(name.color.key[i,5])){
				expr.strain.locale <- which(names(centered.expr[[1]]) == name.color.key[i,4])
				if(length(expr.strain.locale) > 0){
					strain.expr <- unlist(lapply(centered.expr, function(x) x[expr.strain.locale]))[gene.order1]
				
					do.strain.locale <- which(names(coef.list[[1]]) == name.color.key[i,5])
					strain.coef <- unlist(lapply(coef.list, function(x) x[do.strain.locale]))[gene.order2]
					model <- lm(strain.expr~strain.coef)
					r2 <- summary(model)$adj.r.squared
					plot(strain.coef, strain.expr, xlab = "eQTL Coefficient", ylab = "log(Strain Expression) - log(Mean Expression)", main = paste(col.table[i,2], "\nr2 =", signif(r2, 3)), col = name.color.key[i,3])
					abline(model)			
					}
				}
			}
		dev.off()
		
		}

	
	
}