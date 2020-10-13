# performs pairwise comparison of two groups

doEdgeRAnalysis= function(count.table, groups){
	
	y=DGEList(counts=count.table, group=groups)
	
	#filter out low-expressed genes as in the EdgeR user's guide
	good.counts <- rowSums(cpm(y)>1) >= 2
	y <- y[good.counts,,keep.lib.sizes=FALSE]
	y <- calcNormFactors(y) #calculate normalization factors to account for library size
	design <- model.matrix(~groups)
	
	y <- estimateDisp(y,design)
	fit <- glmQLFit(y,design)
	qlf <- glmQLFTest(fit,coef=2)
	
	return(qlf)
	}
	
