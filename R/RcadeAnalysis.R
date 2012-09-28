##Two forms of analysis

##1) Provide args corresponding to necessary resources

##2) TODO Finish "incomplete" analysis - update(Rcade)

##this is 1).

RcadeAnalysis <- function(DE, ChIPannoZones, annoZoneGeneidName, ChIPtargets, ChIPfileDir, cl, DElookup, DE.prior = NULL, prior.mode="assumeIndependent", prior=NULL, ...)
{
	##FIXME allow altering of e.g. getLibsizesArgs
	##TODO ability to dump partial output upon error?
	##TODO use paraList slot - put parameters into slot
	##TODO verbose parsing of ChIPtargets

	##FIXME before we do *anything*, check arguments are valid
	##FIXME check that there is sufficient overlap between geneIDs to proceed
	##ChIPannoZones/annoZoneGeneidName
	##ChIPtargets
	checkTargets(ChIPtargets)
	##ChIPfileDir
	##cl
	if(missing(cl))
	{
		cl <- NULL
	} else {
		if(!is.null(cl) && !"cluster" %in% class(cl)){stop("cl argument must be of 'cluster' class.")}
	}
	##DElookup
	DElookup <- sanitizeDElookup(DE, DElookup)
	##DE
	##DE.prior
	##
	##FIXME NAs? Check column types are correct?

	##---------------------------------------------------------------

	##create new Rcade object, transfer things across
	output <- new("Rcade")
	output@ChIP[[1]]@annoZones <- ChIPannoZones
	output@ChIPtargets <- ChIPtargets
	output@ChIPfileDir <- ChIPfileDir
	output@DE <- DE
	if(!is.null(DE.prior)) output@DE.prior <- DE.prior

	##At this point, we could replace everything with e.g. update(output)??
	##That might be difficult for partial output though.

	##2 countReads -------------------------------
	output@ChIP[[1]]@counts <- countReads(ChIPannoZones, ChIPtargets, ChIPfileDir)

	##3 baySeq -----------------------------------
	output@ChIP[[1]]@summary <- diffCountsBaySeq(output@ChIP[[1]]@counts, ChIPtargets, ChIPannoZones, cl = cl, getLibsizesArgs = list(estimationType = "quantile", quantile = 0.75), getPriors.NBArgs = NULL, getLikelihoods.NBArgs = NULL)
	##FIXME get libsizes from targets
	##FIXME transfer Args
	
	##4 constructRcadeTable ----------------------	
	output@Rcade <- constructRcadeTable(DE, DElookup, output@ChIP[[1]]@summary, ChIPannoZones, annoZoneGeneidName, DE.prior, ChIP.prior=NULL, prior.mode=prior.mode, prior)

	output
}
