##I wish for this class to be "invisible" to most user, hence no accessors
##This may change in future releases
##Please use "Rcade" class accessors/methods

setMethod("show", "RcadeTrack",
	function(object)
		{
			describeGR(object@annoZones, "annoZones")
			cat("shift:"); cat(object@shift, "\n")
			describeMx(object@counts, "counts")
			describeMx(object@summary, "summary")
		}
)
