##Plot methods

##Plot 1: MM plot FIXME
setGeneric("plotMM",
	function(x, DE.abs=FALSE, ...) standardGeneric("plotMM")
)

##Plot 1b: BB plot
setGeneric("plotBB",
	function(x, ...) standardGeneric("plotBB")
)

##Plot 2: 3D B-B-B plot FIXME
setGeneric("plotBBB",
	function(x, ...) standardGeneric("plotBBB")
)


##Plot 3: PCA FIXME
setGeneric("plotPCA",
	function(x, ...) standardGeneric("plotPCA")
)

##Plot 4: View reads in vicinity of gene of interest
##TODO


##Output method

setGeneric("exportRcade",
	function(x, directory="RcadeOutput", cutoffMode="top", cutoffArg = 1000, justGeneID=FALSE, removeDuplicates="beforeCutoff") standardGeneric("exportRcade")
)

##Accessors

setGeneric("getDE", function(x, what="summary") standardGeneric("getDE"))
setGeneric("getChIP", function(x, what="summary") standardGeneric("getChIP"))
setGeneric("getRcade", function(x) standardGeneric("getRcade"))

##Replacers

setGeneric("getDE<-", function(x, what="summary", value) standardGeneric("getDE<-")) 
setGeneric("getChIP<-", function(x, what="summary", value) standardGeneric("getChIP<-")) 
setGeneric("getRcade<-", function(x, value) standardGeneric("getRcade<-")) 
