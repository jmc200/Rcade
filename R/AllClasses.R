## Class for storing a "track" (i.e. data for a ChIP-seq epitope)
setClass("RcadeTrack",
		representation(
			annoZones = "GRanges",
			shift = "integer",
			prior = "numeric",
			counts = "matrix",
			summary = "data.frame"
		),

		prototype(
			shift = 0L
		)
)

## 
setClass("Rcade",
		representation(
			Rcade = "data.frame",
			DE = "data.frame",
			DE.prior = "numeric",
			ChIP = "list", ##N.B. List of RcadeTracks, used to allow compatibility with multiple epitopes. Names correspond to epitope name.
			#ChIPEpitope = "character", - future compatibility
			ChIPtargets = "data.frame",
			ChIPfileDir = "character",
			history = "character",
			paraList = "list"
		),

		prototype(
			ChIP = list(new("RcadeTrack"))
		),

		validity = function(object)
			{
				##require all things in object@ChIP to be RcadeTracks
				if(any(sapply(object@ChIP, class) != "RcadeTrack"))
					return("ChIP slot must be a list of RcadeTrack objects.")
				##TODO require object@ChIP to have epitope names?
				##TODO check references between tables are intact
				TRUE
			}
)
