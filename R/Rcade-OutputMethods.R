##two functions required for FNR + FPR minimization
cumsumbackwards <- function(x)
{
	x <- rev(x)
	x <- cumsum(x)
	x <- rev(x)
}

minFDR <- function(PP, k)
{
	PP <- sort(PP, decreasing = TRUE)
	##k.ind <- 1:n ##how many things should we call?

	FNR <- cumsumbackwards(PP)/(length(PP):1)
	FPR <- cumsum(1-PP)/(1:length(PP))

	FNR <- c(FNR, 0) #call everything, no false negatives
	FPR <- c(0, FPR) #call nothing, no false positives

	#plot(FPR); points(FNR, col = "blue")

	which.min(FNR + k*FPR) - 1 ##number of things that should be called 
}

##useful function for sorting output
outputWriteFunction <- function(x, directory, name, cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by)
{
	##sort, if requested
	if(!is.null(sort.by))
	{
		if(sort.by != "B.ChIP.DE") ##This cheat avoids sorting twice.
		{
			x <- x[order(x[,sort.by], decreasing=TRUE),]
		}
	}

	##Remove duplicates?
	if(removeDuplicates == "beforeCutoff")
	{
		x <- x[!duplicated(x$geneID),]
	}

	##FIXME cutoffMode
	if(tolower(cutoffMode) == "all")
	{
		##no thresholding
	} else if(tolower(cutoffMode) == "top") {
		x <- head(x, cutoffArg)
	} else if(tolower(cutoffMode) == "b") {
		x <- x[x[,sort.by] > cutoffArg,]
	} else if(tolower(cutoffMode) == "fdr") {
		sel <- minFDR(x, cutoffArg)
		x <- head(x, sel)
	} else {
		stop("cutoffMode argument must be 'all', 'top', 'B', or 'FDR'.")
	}

	##Remove duplicates?
	if(removeDuplicates == "afterCutoff")
	{
		x <- x[!duplicated(x$geneID),]
	}

	##output to file
	if(justGeneID)
	{
		write.csv(x$geneID, file = file.path(directory, name), row.names = F)
	} else {
		write.csv(x, file = file.path(directory, name), row.names = F)
	}
}

setMethod("exportRcade", c(x = "Rcade"),
	function(x, directory, cutoffMode, cutoffArg, justGeneID, removeDuplicates)
	{
		#check for valid cutoffMode and removeDuplicates FIXME
		if(!cutoffMode %in% c("all", "top", "B", "FDR"))
			{stop("cutoffMode argument must be one of 'all', 'top', 'B' or 'FDR'.")}
		if(!removeDuplicates %in% c("beforeCutoff", "afterCutoff", "none"))
			{stop("cutoffMode argument must be one of 'beforeCutoff', 'afterCutoff' or 'none'.")}

		temp <- x@Rcade
		temp <- temp[order(temp$B.ChIP.DE, decreasing=TRUE),]

		##check existence of path?
		dir.create(file.path(directory), showWarnings = FALSE)

		##DIRECTLY UPREGULATED
		##only things with DA's M > 0
		##separate files for DE > 0 and DE < 0
		##NB: already sorted by B.ChIP.DE
		sel.up <- temp$M.ChIP > 0 & temp$logfc.DE > 0
		sel.down <- temp$M.ChIP > 0 & temp$logfc.DE < 0
		sel.up <- ifelse(is.na(sel.up), TRUE, sel.up)
		sel.down <- ifelse(is.na(sel.down), TRUE, sel.down)

		outputWriteFunction(temp, directory, "DEandChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE")
		outputWriteFunction(temp[sel.up,], directory, "UpChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE")
		outputWriteFunction(temp[sel.down,], directory, "DownChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE")

		##INDIRECTLY REGULATED
		##as before but now require M < 0 from DA.
		sel.up.np <- temp$M.ChIP < 0 & temp$logfc.DE > 0
		sel.down.np <- temp$M.ChIP < 0 & temp$logfc.DE < 0
		sel.up.np <- ifelse(is.na(sel.up.np), TRUE, sel.up.np)
		sel.down.np <- ifelse(is.na(sel.down.np), TRUE, sel.down.np)

		outputWriteFunction(temp[sel.up.np,], directory, "UpNoChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE")
		outputWriteFunction(temp[sel.down.np,], directory, "DownNoChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE")

		##UP/DOWNREGULATED
		##things with highest B.DE
		##separate files for DE > 0 and DE < 0
		sel.up <- temp$logfc.DE > 0
		sel.down <- temp$logfc.DE < 0

		outputWriteFunction(temp[sel.up,], directory, "Up.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.DE")
		outputWriteFunction(temp[sel.down,], directory, "Down.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.DE")

		##~ "POISED"
		##things with binding but no DE
		sel <- temp$M.ChIP > 0
		outputWriteFunction(temp[sel,], directory, "ChIPonly.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.only")

		##ChIP
		##things with binding
		outputWriteFunction(temp[sel,], directory, "ChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP")

		##WORST
		##things with absolutely nothing
		outputWriteFunction(temp, directory, "Nothing.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.nothing")
	
	}
)

