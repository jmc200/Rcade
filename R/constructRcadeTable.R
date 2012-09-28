constructRcadeTable <- function(DE, DElookup, chip, annoZone, annoZoneGeneidName, DE.prior=NULL, ChIP.prior=NULL, prior.mode, prior=NULL)
{
	##FIXME What happens when we have incomplete data? e.g. blank geneIDs
	##FIXME allow choice of prior; error on missing information

	##FIXME check prior
	if (prior.mode == "keepChIP") {
		if(any(is.null(c(DE.prior, prior))))
		{
			stop("For prior.mode='keepChIP', must specify both 'DE.prior' and 'prior' arguments.")
		}
	}

	##start with ChIP
	##append annotation information to ChIP
	chip$geneID <- elementMetadata(annoZone)[,1] ##FIXME check for breakability
	##Optional - mass annotation TODO (can remove this?)
	#ord <- as.numeric(rownames(chip))
	#chip <- cbind(chip, as.data.frame(TSS)[ord,])

	##sort table
	chip <- chip[order(chip$B, decreasing=T),]

	##DE--------------------------------------------------------------------------
	##sanitize
	DElookup <- sanitizeDElookup(DE, DElookup)

	##Extract relevant DE cols
	DElookupV <- unlist(DElookup)
	DE.sel <- DE[,DElookupV]
	colnames(DE.sel) <- names(DElookupV) ##use canonical colnames

	##combination of DE and DA information ---------------------------------------

	##Merge DA and DE information
	colnames(DE.sel) <- paste(colnames(DE.sel), ".DE", sep = "") ##force suffix
	colnames(chip) <- paste(colnames(chip), ".ChIP", sep = "") ##force suffix
	rcade <- merge(DE.sel, chip, by.x = "geneID.DE", by.y = "geneID.ChIP")
	colnames(rcade)[1] <- "geneID"

	##Joint analysis
	##calculate relevant probabilities
	rcade$p.DE <- logit(rcade$B.DE)
	rcade$p.ChIP <- logit(rcade$B.ChIP)

	##FIXME considerations re: prior
	if(prior.mode == "assumeIndependent")
	{
		rcade$B.nothing <- with(rcade, arclogit((1 - p.DE) * (1 - p.ChIP)))
		rcade$B.DE.only <- with(rcade, arclogit(p.DE * (1 - p.ChIP)))
		rcade$B.ChIP.only <- with(rcade, arclogit((1 - p.DE) * p.ChIP))
		rcade$B.ChIP.DE <- with(rcade, arclogit(p.DE * p.ChIP))
	} else if (prior.mode == "keepChIP") {
		names(prior) <- c("D|C", "D|-C") ##this is mostly for my benefit, really
		probs <- rep(0, 4)
		names(probs) <- c("-C-D", "-CD", "C-D", "CD")
		probs <- with(rcade, list(
			"CD"=p.DE*p.ChIP*prior[1]/DE.prior,
			"-CD"=p.DE*(1-p.ChIP)*prior[2]/DE.prior,
			"C-D"=(1-p.DE)*p.ChIP*(1-prior[1])/(1-DE.prior),
			"-C-D"=(1-p.DE)*(1-p.ChIP)*(1-prior[2])/(1-DE.prior)
		))
		probs <- probs/sum(probs)
		rcade$B.nothing <- arclogit(probs[["-C-D"]])
		rcade$B.DE.only <- arclogit(probs[["-CD"]])
		rcade$B.ChIP.only <- arclogit(probs[["C-D"]])
		rcade$B.ChIP.DE <- arclogit(probs[["CD"]])
	}

	rcade <- rcade[order(rcade$B.ChIP.DE, decreasing = T),]

	##Annotation
	#rcade$Genes.Symbol <- convertENSGtoOGS(o$ENSG, mart)

	rcade
}
