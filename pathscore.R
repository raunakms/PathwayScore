### FUNCTON TO COMPUTE PATHWAY Z-SCORE --------------------------------------------------------------
# Usage:
# get.pathway.zcore(file.expr, file.pathway)
#
# Description:
# file.expr = path of file containing gene-expression matrix with Row as Genes and Columns as Samples.
# file.pathway = path of file containing pathway-gene relationship
#	- each line should contain information about a unique pathway
# 	- first-column must contain name of the pathway
#	- second-column must contain genes mapped to the corresponding pathway 
#		multiple genes must be separated by comma symbol (,)
### --------------------------------------------------------------------------------------------------
get.pathway.zcore <- function(file.expr, file.pathway){
	# LOAD LIBRARIES -----
	require("stringr")

	# LOAD Expression ----
	expr <- read.delim(file.expr, header=T, row.names=1, stringsAsFactors=F)

	# LOAD Pathway ------
	dat.path <- read.delim(file.pathway, header=T, stringsAsFactors=F)

	# Get Zscores -------
	pathz <- matrix(0, nrow=nrow(dat.path), ncol=4, dimnames=list(dat.path[,1], colnames(expr)[-1]))
	for(i in 1:nrow(dat.path)){
		pathway <- dat.path[i,1]
		genes <- str_split(dat.path[i,2], ",")[[1]]

		expr.path <- subset(expr, rownames(expr) %in% genes)[,-1]
		exprz <- getZscore(expr.path)
		pathz[i,] <- apply(exprz, 2, mean)
	}
	return(pathz)
}


### FUNCTION TO CALCULATE Z-SCORE ----------------------------------------
# Let, M be the gene-expression matrix 
# with Row as Genes and Columns as Samples.
#
# Zscore = (x_ij - m_i) / sd_i
# Here,
# i = gene
# j = sample
# x_ij = expression value of gene i in sample j
# m_i = mean of expression values in for gene i across all j
# sd_i = standard-deviation of expression values in for gene i across all j
#
# Usage:
# getZscore(dat)
#
# Description:
# dat = gene-expression matrix with Row as Genes and Columns as Samples.
### ------------------------------------------------------------------------
getZscore <- function(dat){
	z.dat <- matrix(nrow=nrow(dat), ncol=ncol(dat), dimnames=list(rownames(dat),colnames(dat)))

	dat.mean <- apply(dat, 1, mean)
	dat.stdev <- apply(dat, 1, sd)

	for(i in 1:nrow(dat)){
		x <- as.numeric(dat[i,])
		z.dat[i,] <- (x - dat.mean[i])/dat.stdev[i]
	}
	z.dat <- as.data.frame(z.dat)
	return(z.dat)
}
