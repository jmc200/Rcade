##Plot 1: MM plot
setMethod("plotMM", c(x = "Rcade"),
	function(x, DE.abs, xlab="DE Log Ratio", ylab="ChIP Log Ratio", pch=19, xlim=NULL, ylim=NULL, zlim=NULL, main="ChIP against Expression", col.scale=rainbow(10000)[1000:10000], col.z.inf="black", legend=TRUE, legend.x=NULL, ...)
	{
	  if(legend)
    {
      require(plotrix) ##need color.legend()
	  }
    
    
		DE <- x@Rcade$logfc.DE
		if(DE.abs)
		{
			DE <- abs(DE)
		}

		if(is.null(xlim))
		{
			temp <- DE[is.finite(DE)]
			xlim <- range(temp)
      if(legend) {xlim[2] <- xlim[2]*8/7.3 + 0.6} ##make room for legend
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

    col.which <- floor(length(col.scale)*heights)
    col.which <- ifelse(col.which < length(col.scale), col.which + 1, col.which) ##i.e. add one, UNLESS heights == 1 (zero-index issue)
    
    cols <-  col.scale[col.which]
    cols <- ifelse(is.na(cols), col.z.inf, cols) ##replace Infs (or any other problematic "heights" value) s.t. col=col.z.inf
    
		#sel <- x@Rcade$M.ChIP > 0

		plot(DE, x@Rcade$M.ChIP, xlab=xlab, ylab=ylab, pch=pch, col=cols, xlim=xlim, ylim=ylim, main=main, ...)
    
		if(legend)
		{	#construct legend
		  legend <- seq(from=zlim[1], to=zlim[2], length.out=6)
		  legend <- round(legend, digits=1)
		  legend <- paste("B = ", legend, sep = "")
		  
		  if(is.null(legend.x))
        {
          legend.x <- xlim[1] + (xlim[2]-xlim[1])*c(7.3,7.4)/8 ##TODO improve? How wide is the text?
		    }
		  
		  ##reduce resolution of colour scale to stop white flecks appearing
		  scale.sel <- floor(seq(from=1, to=length(col.scale), length.out=500))
		  col.scale.sel <- col.scale[scale.sel]
      
		  color.legend(legend.x[1], ylim[1], legend.x[2], 0.15*ylim[1] + 0.85*ylim[2],
		               legend=legend, rect.col=col.scale.sel, cex=1, align="rb", gradient="y")
		  color.legend(legend.x[1], 0.1*ylim[1] + 0.9*ylim[2], legend.x[2], ylim[2],
		               legend="B = Inf", rect.col=col.z.inf, cex=1, align="rb", gradient="y")
		}
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





