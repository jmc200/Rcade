sanitizeDElookup <- function(DE, DElookup)
{
	##The outcome of this function should be a "nice" DElookup object:
	##a list, elements corresponding to targets columns
	##and names correspond to Rcade's names,

	##1) Coerce into required format
	if(missing(DElookup)) {DElookup = list()} ##we'll be OK as long as the col names are obvious
	DE <- as.list(DE)
	if(is.null(names(DElookup)))
	{
		names(DElookup) = DElookup ##if names blank, add them
	} else { 
		sel <- names(DElookup) %in% c("", NA) ##find unnamed fields...
		names(DElookup)[sel] = DElookup[sel] ##...replace with obvious name choices
	}

	names(DElookup) <- tolower(names(DElookup)) ##TODO alternative policy on case?

	##2) Stop if there are duplicates in DElookup
	sel <- duplicated(tolower(names(DElookup)))
	if(any(sel))
	{
		stop("Duplicates present in DElookup:", paste(names(DElookup)[sel], collapse = ", "), ". Please review DElookup object." )
	}

	##3) Check for required entries - cannot be missing or duplicated 
	required <- c("geneid", "logfc", "b") ##<<<--------------------------------REQUIRED

	##which are not present in DElookup?
	missing <- is.na(pmatch(required, tolower(names(DElookup))))
	##directly present in DE?
	if(any(missing))
	{
		sel <- pmatch(required[missing], tolower(names(DE))) ##only look for non-missing fields
		missingDE <- is.na(sel)

		##if any fields are missing here, we are dead in the water
		if(any(missingDE)) {stop("DElookup missing required field(s): ", paste(required[missing][missingDE], collapse = ", "))}

		##insert missing fields
		temp <- as.list(names(DE)[sel])
		names(temp) <- temp
		DElookup <- c(DElookup, temp)
	}

	##guess missing entries of lookup table? FIXME
	##...
	##USEFUL: probeID, t value, other annotation information...

	#DElookup <- data.frame(ENSG = ENSG.col, coef = coef.col, t = t.col, p.adj = p.adj.col, B = B.col, coef.vec = coef.vec.col, t.vec = t.vec.col, p.adj.vec = p.adj.vec.col, B.vec = B.vec.col, ProbeID = ProbeID.col)
	names(DElookup) <- tolower(names(DElookup)) ##TODO alternative policy on case?
	##return "b" -> "B" and "geneid" -> "geneID", for consistancy with ChIP data
	names(DElookup)[names(DElookup) == "b"] <- "B"
	names(DElookup)[names(DElookup) == "geneid"] <- "geneID"

	DElookup
}
