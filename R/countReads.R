##get the levels in the Bam file
getBamChrs <- function(file)
{
	require(Rsamtools)
	header <- scanBamHeader(file)
	header <- lapply(header[[1]]$text, function(x) {x})
	header <- header[names(header) == "@SQ"] ##get sequence lines only
	header <- sapply(header, function(x){x[grep("SN:", x)]})##get sequence names only
	header <- substr(header, 4, nchar(header))
	header
}

##----------------------------------------------------------
##Construct count table

countReads <- function(annoZone, targets, fileDir = NULL, dontCheckTargets=FALSE) ##targets must have columns: fileID, sampleID, factor, filepath
{
	##TODO support for things that aren't bams

	##interpret "targets" argument & various other checks
	if(!dontCheckTargets)
	{
		targets <- checkTargets(targets)
	}

	#check existence of all files before continuing
	files <- file.path(fileDir, as.character(targets$filepath))
	sel <- file.exists(files)

	if(any(!sel))
	{
		stop("Cannot find files: ", paste(files[!sel], collapse = ", "))
	}

	##check for index files
	if("index" %in% colnames(targets))
	{
		index <- file.path(fileDir, as.character(targets$index))
	} else {
		index <- files
	}
	##if we can't find the index files, then don't use any
	index.real <- paste(index, ".bai", sep="")
	sel <- file.exists(index.real)
	if(any(!sel))
	{
		message("Cannot find index files: ", paste(index.real, collapse = ", "))
		message("Therefore, I will not use any index files") ##FIXME Just for the missing ones?
		index <- NULL
	}

	##BAM file readin parameters
	what <- c("rname", "strand", "pos", "qwidth")
	p <- ScanBamParam(what = what)

	##prepare output matrix
	output <- matrix(0, nrow = length(annoZone), ncol = nrow(targets), dimnames = list(rownames(annoZone), targets$sampleid))

	##-------------------------------------------

	##Preliminary checks - any unfixable errors?
	annoChrs <- seqlevels(annoZone)
	annoChrsChomp <- gsub("chr", "", annoChrs, ignore.case = TRUE)

	chompflag <- rep(FALSE, ncol(output))
	for(i in 1:ncol(output)) #TODO rewrite without for loop
	{
		file <- as.character(targets$filepath[i])

		bamChrs <- getBamChrs(files[i])

		##Check for chromosome name mismatch (usually just "chr" presence/absence)
		message("File ", file, ": Found ", length(bamChrs), " chromosomes, of which ", sum(bamChrs %in% annoChrs), " are present in the annotation.")
		if(all(!bamChrs %in% annoChrs))
		{
			##if we were to remove all instances of "chr", does this fix things?
			bamChrsChomp <- gsub("chr", "", bamChrs, ignore.case = TRUE)
			if(all(!bamChrsChomp %in% annoChrsChomp))
			{
				stop("File ", file, ": Could not match chromosomes with annotation. Chromosomes:", bamChrs,)
			} else {
				message("File ", file, ": Removed all instances of 'chr' - now ", sum(bamChrsChomp %in% annoChrsChomp), " chromosomes are present in the annotation.")
				#message("File ", file, ": ", !sum(bamChrs %in% annoChrs), " chromosomes still not present in annotation.")
				chompflag[i] <- TRUE
				#print(bamChrs) -- diagnostic
			}
		}
	}

	##-------------------------------------------

	##Read in each sample...
	for(i in 1:ncol(output)) #TODO rewrite without for loop
	{
		file <- as.character(targets$filepath[i])
		message("File ", file, ": Counting reads.")

		#read in file
		if(is.null(index))
		{
			x <- scanBam(file=files[i], param=p)[[1]]
		} else {
			x <- scanBam(file=files[i], index=index[i], param=p)[[1]]
		}

		##TODO Alternative: only collect reads of interest. Currently slower, in theory should be faster? Certainly better wrt memory.
#		p.derp <- ScanBamParam(what = what, which = annoZone)

#		x <- scanBam(file = files[i], index = files[i], param = p.derp)

		##FIXME remove NAs (probably due to QC leaving malformed reads?)
		if(any(is.na(x$rname)) | any(is.na(x$pos)))
		{
			message("File ", file, ": File contains NAs. (Rcade has removed them.)")
			warning("File ", file, ": File contains NAs. (Rcade has removed them.)")
			sel <- !is.na(x$rname) & !is.na(x$pos)
			x <- lapply(x, function(a) {a[sel]})
		}

		##construct GRanges
		x <- with(x,
			GRanges(
				seqnames = as(rname, "Rle"),
				ranges = IRanges(start = pos, width = qwidth),
				strand = as(strand, "Rle")
			)
		)

		##remove chr if appropriate
		if(chompflag[i])
		{
			seqlevels(annoZone) <- annoChrsChomp ##NB case throughout function
			seqlevels(x) <- gsub("chr", "", seqlevels(x), ignore.case = TRUE)
		} else {
			seqlevels(annoZone) <- annoChrs
		}

#		##make sure both TSS and x are on the same spaces
#		if(any(!names(x) %in% names(annoZone)))
#		{
#			sel <- names(x)[names(x) %in% names(annoZone)]
#			x <- x[sel]
#		}

		ChIPshift = 0 ##FIXME obtain this from targets

		##reduce to 5'
		sel <- strand(x) == "+"
		start <- IRanges::ifelse(sel, start(x) + ChIPshift, end(x) - ChIPshift)
		start(x) <- start
		width(x) <- 1

		##kill missing spaces(worth doing?) TODO
		
		##get the bin counts
		counts <- countOverlaps(annoZone, x)

		##update output
		output[,i] <- counts

		rm(x); gc()
	}

	colnames(output) <- targets$sampleid

	output
}
