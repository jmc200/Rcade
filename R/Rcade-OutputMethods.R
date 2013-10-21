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
outputWriteFunction <- function(x, directory, name, cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by, x.sel.after=NULL)
{
	##append x.sel.after to x, (resolved later)
	if(!is.null(x.sel.after))
	{
		if(length(x.sel.after) != nrow(x) | !is.logical(x.sel.after))
		{
			stop("x.sel.after must be a TRUE/FALSE vector with length equal to nrow(x)")
		}
		x <- cbind(x, PUTNJZYOGKHMDWQLEAFXSBCVRI=x.sel.after) ##random column name, must not == user column name
	}

	##sort, if requested
	if(!is.null(sort.by))
	{
		if(sort.by != "B.ChIP.DE") ##No need to sort. Already sorted by this value
		{
			x.order <- order(x[,sort.by], decreasing=TRUE)
			x <- x[x.order,]
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
		sel <- minFDR(expit(x[,sort.by]), cutoffArg)
		x <- head(x, sel)
	} else {
		stop("cutoffMode argument must be 'all', 'top', 'B', or 'FDR'.")
	}

	##Remove duplicates?
	if(removeDuplicates == "afterCutoff")
	{
		x <- x[!duplicated(x$geneID),]
	}

	##take only the rows with x.sel.after == TRUE (if appropriate)
	if(!is.null(x.sel.after))
	{
		x <- x[x[,"PUTNJZYOGKHMDWQLEAFXSBCVRI"], ]
		x <- x[,colnames(x) != "PUTNJZYOGKHMDWQLEAFXSBCVRI"]
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
		sel <- temp$M.ChIP > 0
		sel.up <- sel & temp$logfc.DE > 0
		sel.down <- sel & temp$logfc.DE < 0
		sel <- ifelse(is.na(sel), TRUE, sel)
		sel.up <- ifelse(is.na(sel.up), TRUE, sel.up)
		sel.down <- ifelse(is.na(sel.down), TRUE, sel.down)

		outputWriteFunction(temp, directory, "DEandChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE", x.sel.after=sel)
		outputWriteFunction(temp, directory, "UpChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE", x.sel.after=sel.up)
		outputWriteFunction(temp, directory, "DownChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE", x.sel.after=sel.down)

		##INDIRECTLY REGULATED
		##as before but now require M < 0 from DA.
		sel.up.np <- temp$M.ChIP < 0 & temp$logfc.DE > 0
		sel.down.np <- temp$M.ChIP < 0 & temp$logfc.DE < 0
		sel.up.np <- ifelse(is.na(sel.up.np), TRUE, sel.up.np)
		sel.down.np <- ifelse(is.na(sel.down.np), TRUE, sel.down.np)

		outputWriteFunction(temp, directory, "UpNoChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE", x.sel.after=sel.up.np)
		outputWriteFunction(temp, directory, "DownNoChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.DE", x.sel.after=sel.down.np)

		##UP/DOWNREGULATED
		##things with highest B.DE
		##separate files for DE > 0 and DE < 0
		sel.up <- temp$logfc.DE > 0
		sel.down <- temp$logfc.DE < 0

		outputWriteFunction(temp, directory, "Up.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.DE", x.sel.after=sel.up)
		outputWriteFunction(temp, directory, "Down.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.DE", x.sel.after=sel.down)

		##~ "POISED"
		##things with binding but no DE
		sel <- temp$M.ChIP > 0
		outputWriteFunction(temp, directory, "ChIPonly.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP.only", x.sel.after=sel)

		##ChIP
		##things with binding
		outputWriteFunction(temp, directory, "ChIP.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.ChIP", x.sel.after=sel)

		##WORST
		##things with absolutely nothing
		outputWriteFunction(temp, directory, "Nothing.csv", cutoffMode, cutoffArg, justGeneID, removeDuplicates, sort.by="B.nothing")
	
	}
)

