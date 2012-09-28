##convenience functions

describeMx <- function(x, nm)
{
	if(is.null(x)) {cat(nm, ": NULL", sep = "")} 
	else cat(nm, ": ", nrow(x), " rows, ", ncol(x), " columns.\n", sep = "")
}

describeVec <- function(x, nm)
{
	if(is.null(x)) {cat(nm, ": NULL", sep = "")} 
	else cat(nm, ": length ", length(x), ".\n", sep = "")
}

describeGR <- function(x, nm)
{
	if(is.null(x)) {cat(nm, ": NULL", sep = "")} 
	else
	{
		lx <- length(x)
		nc <- ncol(elementMetadata(x))
		cat(nm, ": ", class(x), " with ", lx, " ", ifelse(lx == 1L, "range", 
			"ranges"), " and ", nc, " elementMetadata ", ifelse(nc == 
			1L, "col", "cols"), ".\n", sep = "")
	}
}

##show

setMethod("show", "Rcade", 
	function(object)
		{
			describeMx(object@Rcade, "Rcade")
			describeMx(object@DE, "DE")

			##ChIP
			cat("ChIP: ", length(object@ChIP), " track(s):\n", sep = "")
			#list tracks
			if(length(object@ChIP) > 0)
			{
				for(i in 1:length(object@ChIP))
				{
					show(object@ChIP)
				}
			}

			cat("Metadata:\n")
			describeMx(object@ChIPtargets, "ChIPtargets")
			cat("fileDir:", object@ChIPfileDir, "\n", sep = "")
			#describeVec(object@history, "history") For future compatability
			#describeVec(object@paraList, "paraList") For future compatability
		}
)

##accessors
setMethod("getDE", c(x="Rcade"), 
	function(x, what)
		{
			if(what == "summary")
			{
				return(x@DE)
			} else if (what == "prior") {
				return(x@DE.prior)
			} else {
				stop('"what" must be "summary" or "prior".')
			}
		}
)

setMethod("getChIP", c(x="Rcade"), 
	function(x, what)
		{
			if(what == "targets")
			{
				return(x@ChIPtargets)
			} else {
				return(slot(x@ChIP[[1]], what))
			}
		}
)

setMethod("getRcade", c(x="Rcade"), 
	function(x) {x@Rcade}
)

##replacers
setReplaceMethod("getDE", c(x="Rcade"), 
	function(x, what, value)
		{
			if(what == "summary")
			{
				x@DE <- value
			} else if (what == "prior") {
				x@DE.prior <- value
			} else {
				stop('"what" must be "summary" or "prior".')
			}
		}
)

setReplaceMethod("getChIP", c(x="Rcade"), 
	function(x, what, value)
		{
			if(what == "targets")
			{
				x@ChIPtargets <- value
			} else {
				slot(x@ChIP[[1]], what) <- value
			}
			x
		}
)

setReplaceMethod("getRcade", c(x="Rcade"), 
	function(x, value)
		{
			x@Rcade <- value
			x
		}
)
