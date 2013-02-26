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
	rcade$p.DE <- expit(rcade$B.DE)
	rcade$p.ChIP <- expit(rcade$B.ChIP)

	##FIXME considerations re: prior
	if(prior.mode == "assumeIndependent")
	{
		rcade$B.nothing <- with(rcade, logit((1 - p.DE) * (1 - p.ChIP)))
		rcade$B.DE.only <- with(rcade, logit(p.DE * (1 - p.ChIP)))
		rcade$B.ChIP.only <- with(rcade, logit((1 - p.DE) * p.ChIP))
		rcade$B.ChIP.DE <- with(rcade, logit(p.DE * p.ChIP))
	} else if (prior.mode == "keepChIP") {
		names(prior) <- c("D|C", "D|-C") ##this is mostly for my benefit, really
		probs <- rep(0, 4)
		names(probs) <- c("nCnD", "nCD", "CnD", "CD") #'-' not allowed in data.frame names? replaced with n
		probs <- with(rcade, data.frame(
			"CD"=p.DE*p.ChIP*prior[1]/DE.prior,
			"nCD"=p.DE*(1-p.ChIP)*prior[2]/DE.prior,
			"CnD"=(1-p.DE)*p.ChIP*(1-prior[1])/(1-DE.prior),
			"nCnD"=(1-p.DE)*(1-p.ChIP)*(1-prior[2])/(1-DE.prior)
		))

		probs <- probs/rowSums(probs)


		rcade$B.nothing <- logit(probs[["nCnD"]])
		rcade$B.DE.only <- logit(probs[["nCD"]])
		rcade$B.ChIP.only <- logit(probs[["CnD"]])
		rcade$B.ChIP.DE <- logit(probs[["CD"]])
	}

	rcade <- rcade[order(rcade$B.ChIP.DE, decreasing = T),]

	##Annotation
	#rcade$Genes.Symbol <- convertENSGtoOGS(o$ENSG, mart)

	rcade
}
