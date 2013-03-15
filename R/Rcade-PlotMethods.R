##Plot 1: MM plot
##TODO allow |logfc|
setMethod("plotMM", c(x = "Rcade"),
	function(x, DE.abs, xlab="DE Log Ratio", ylab="ChIP Log Ratio", pch=19, xlim=NULL, ylim=NULL, zlim=NULL, main="ChIP against Expression", col.scale=rainbow(10000), ...)
	{
		DE <- x@Rcade$logfc.DE
		if(DE.abs)
		{
			DE <- abs(DE)
		}

		if(is.null(xlim))
		{
			temp <- DE[is.finite(DE)]
			xlim <- range(temp)
		}

		if(is.null(ylim))
		{
			temp <- x@Rcade$M.ChIP[is.finite(x@Rcade$M.ChIP)]
			ylim <- range(temp)
			ylim[1] <- 0
		}

		heights <- x@Rcade$B.ChIP.DE
		if (is.null(zlim)) {
			zlim <- range(heights, na.rm=TRUE, finite=TRUE)
		}

        heights <- heights - zlim[1]
        heights <- heights/(zlim[2] - zlim[1])

		#sel <- x@Rcade$M.ChIP > 0

		plot(DE, x@Rcade$M, xlab=xlab, ylab=ylab, pch=pch, col=col.scale[floor(length(col.scale)*heights) + 1], xlim=xlim, ylim=ylim, main=main, ...)
	}
)

##Plot 1b: BB plot
setMethod("plotBB", c(x = "Rcade"),
	function(x, xlab="DE Log-odds", ylab="ChIP Log-odds", pch=19, xlim=NULL, ylim=NULL, main="ChIP against Expression", ...)
	{
		DE <- x@Rcade$B.DE
		ChIP <- x@Rcade$B.ChIP

		if(is.null(xlim))
		{
			temp <- DE[is.finite(DE)]
			xlim <- range(temp)
		}
		if(is.null(ylim))
		{
			temp <- ChIP[is.finite(ChIP)]
			ylim <- range(temp)
		}

		plot(DE, ChIP, xlab=xlab, ylab=ylab, pch=pch, xlim=xlim, ylim=ylim, main=main, ...)
	}
)

##Plot 2: 3D B-B-B plot
setMethod("plotBBB", c(x = "Rcade"),
	function(x, xlab = "B (DE)", ylab = "B (ChIP)", zlab = "B (Combined)", main = "ChIP against Expression", col.scale=rainbow(10000), col, ...)
	{
		require(rgl)

		if(missing(col))
		{
			heights <- x@Rcade$B.ChIP.DE
			zlim <- range(heights, na.rm=TRUE, finite=TRUE)
		    heights <- heights - zlim[1]
		    heights <- heights/(zlim[2] - zlim[1])

			col=col.scale[length(col.scale)*(heights)]
		}

		sel <- is.infinite(x@Rcade$B.ChIP.DE)

		if(any(sel))
		{
			warning("Removed at least one gene with infinite B.ChIP.DE")
			plot3d(x@Rcade$B.DE[!sel], x@Rcade$B.ChIP[!sel], x@Rcade$B.ChIP.DE[!sel], xlab=xlab, ylab=ylab, zlab=zlab, col=col[!sel], ...)
			decorate3d(main = main)
		} else {
			plot3d(x@Rcade$B.DE, x@Rcade$B.ChIP, x@Rcade$B.ChIP.DE, xlab=xlab, ylab=ylab, zlab=zlab, col=col, ...) ##or heat.colors?
			decorate3d(main = main)
		}
	}
)


##Plot 3: PCA

setMethod("plotPCA", c(x = "Rcade"),
	function(x, col = "black", pch = NULL, main = "ChIP PCA", ...)
	{
		cs <- x@ChIP[[1]]@counts

		pc <- prcomp(t(cs)/colSums(cs))

		if(is.null(col))
		{
			##FIXME automagically style colours by sample, etc...
		}

		if(is.null(pch))
		{
			pch = as.character(1:ncol(x@ChIP[[1]]@counts))
		}

		plot(pc$x, main=main, col=col, pch=pch, ...)

		#cols <- c("Red", "Blue", "Pink", "Cyan", "Black")
		#names(cols) <- c("Apo", "Qui", "Sen", "Checkpoint", "Gro")
		#reps <- as.numeric(as.character(targets$Replicate))

		#plot(pc$x, col = cols[as.character(targets$Condition)], pch = ifelse(targets$Factor == "p53", reps, "C"), main = "TSS zone counts - PCA")
		#cols2 <- c("Red", "Blue", "Pink", "Cyan", "Black")
		#names(cols2) <- c("Apoptosis", "Quiescence", "Senescence", "Checkpoint", "Growing")
		#legend(-0.008, 0, legend = names(cols2), col = cols2, pch = "X")
	}
)

##Plot 4: View reads in vicinity of gene of interest TODO





