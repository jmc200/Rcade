##1) must be able to recreate count data from .bam files 

test_countReads <- function() {
	##reference
	data(RcadeSTAT1)
	counts.ref <- getChIP(RcadeSTAT1, what="counts")

	##calculated
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	targets <- read.csv(file.path(dir, "targets.csv"), as.is = TRUE)

	anno <- read.csv(file.path(dir, "anno.csv"))
	anno <- anno[order(anno$chromosome_name),]
	colnames(anno) <- c("ENSG","chr","start","end","str")
	ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")

	counts.cal <- countReads(ChIPannoZones, targets, fileDir = dir)

	##check
	checkEquals(counts.ref, counts.cal)
}

##2) check that annotation zones match between Rcade object and anno.csv

test_RcadeSTAT1 <- function() {
	data(RcadeSTAT1)
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	annoZone.ref <- getChIP(RcadeSTAT1, what="annoZones")

	anno <- read.csv(file.path(dir, "anno.csv"))
	anno <- anno[order(anno$chromosome_name),]
	colnames(anno) <- c("ENSG","chr","start","end","str")

	annoZone.cal <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")

	checkIdentical(annoZone.ref, annoZone.cal)
}

##3) must be able to recreate Rcade table from individual components

test_constructRcadeTable <- function() {
	
	##reference
	data(RcadeSTAT1)
	r.ref <- getRcade(RcadeSTAT1)

	##calculation
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	DE <- getDE(RcadeSTAT1)
	DElookup <- list(GeneID="ENSG", logFC="logFC", B="B", "Genes.Location", "Symbol")
     
	chip <- getChIP(RcadeSTAT1)
	annoZone <- getChIP(RcadeSTAT1, what="annoZones")

	###################################################
  ##prior used in vignette
  ###################################################
	DE.prior = 0.01
	prior.mode = "keepChIP"
	prior = c("D|C" = 0.05, "D|notC" = 0.005)
	cl <- NULL

  r.cal <- constructRcadeTable(DE, DElookup, chip, annoZone, annoZoneGeneidName="ENSG", DE.prior=DE.prior, ChIP.prior=NULL, prior.mode=prior.mode, prior=prior)
  
  #RcadeSTAT1 <- RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName="ENSG",
	#                            ChIPtargets=targets, ChIPfileDir = dir,
	#                            cl=cl, DE.prior=DE.prior, prior.mode=prior.mode, prior=prior,
	#                            DElookup=DElookup)
	#r.cal <- constructRcadeTable(DE, DElookup, chip, annoZone, annoZoneGeneidName="ENSG", DE.prior=NULL, ChIP.prior=NULL, prior.mode="assumeIndependent", prior=NULL)

	cols <- c("B.nothing", "B.DE.only","B.ChIP.only","B.ChIP.DE")

	checkEqualsNumeric(r.ref[,cols], r.cal[,cols])
}

##4) TODO Check baySeq calculation hasn't changed

# test_diffCountsBaySeq <- function() {
#   
#   ##reference
#   data(RcadeSTAT1)
#   r.ref <- ##FIXME get from ... (RcadeSTAT1)
#   
#   ##calcuation
#   woir <- diffCountsBaySeq(RcadeSTAT1@ChIP[[1]]@counts, RcadeSTAT1@ChIPtargets, annoZones = RcadeSTAT1@ChIP[[1]]@annoZones)
# }

##5) TODO Check the prior calculations shown in the vignette


