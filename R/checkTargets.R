checkTargets <- function(targets)
{
	##Check that targets file is appropriate i.e. correct columns, correct entries
	##Cols needed:
	required <- c("fileid", "sampleid", "factor", "filepath")

	targets <- as.data.frame(targets)
	colnames(targets) <- tolower(colnames(targets)) ##force lower case (makes things much easier)
	map <- pmatch(required, tolower(colnames(targets)), nomatch = NA_integer_, duplicates.ok = FALSE)
	names(map) = required

	##FIXME optional columns: shift, scale, index

	##Disallow duplicated col names
	if(any(duplicated(colnames(targets))))
	{
		stop("targets - Duplicated column names not allowed (ignoring case).")
	}

	##------------------------------------------------
	##CHECKS:

	##Disallow if any required column is unmatched
	if(any(is.na(map)))
	{
		stop("targets - could not find match for the following column(s): ", paste(required[is.na(map)], collapse = ", "))
	}

	##check fileIDs are all distinct
	if(any(duplicated(targets[,map["fileid"]]))) {stop("targets - Duplicate fileIDs present.")}

	##Input column present?
	test <- unique(tolower(targets[,map["factor"]]))
	if(!"input" %in% test) {stop("Could not find any 'input' entries in targets$factor.")}
	if(length(test) == 1) {stop("Could not find any non-input entries in targets$factor.")}
	if(length(test) > 2) {stop("Rcade expects 1 non-input entry in targets$factor - found", length(test) - 1, ".")} ##TODO support for multiple factors

	##Samples must not straddle factor levels.
	sel <- tolower(targets[,map["factor"]]) == "input"
	test <- intersect(targets[sel,map["sampleid"]], targets[!sel,map["sampleid"]])
	if(length(test) > 0) {stop(test, ": A sample cannot have multiple factors in 'targets'. ")}

	##------------------------------------------------

	targets
}
