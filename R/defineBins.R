## (STEP 1) Define bins based on some annotation

##example imput
##

#annoTest <- getBM(
#		attributes= c("ensembl_gene_id", "chromosome_name", "transcript_start", "transcript_end", "strand"),
#		mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#	)


##

defineBins <- function(anno, zone, geneID = "ensembl_gene_id", removeDuplicates=TRUE)
{
	##interpret annotation (assume data frame for now? i.e. straight output of biomaRt)
	anno <- as.data.frame(anno)

	##FIXME check that geneID is OK
	if(!geneID %in% colnames(anno))
	{
		stop("geneID not found in colnames(anno).")
	}
	
	##find str, chr, start, end, geneID columns
	##FIXME non-biomaRt input
	chr = "chr"
	start = "start"
	end = "end"
	str = "str"

	##check for "*" problems (needs to work for factors, non-factors)
	if(sum(zone) != 0)	##only a problem if the zone is asymmetric
	{
		temp <- unique(anno[,str]) ##present strand levels...
		sel <- !temp %in% c(1,-1,"+","-") ##unknown strand levels...
		if(any(sel))
		{
			message(paste("Annotation strand information contained the following characters, presumably corresponding to unknown strand:", temp[sel], sep = " "))
			message("I am discarding entries with these strand values.")

			anno <- anno[anno[,str] %in% c(1,-1,"+","-"),]
		}
	}

	##TODO sort data frame? (improves memory usage, but should be user's choice?)

	anno[,start] <- as.integer(ifelse(anno[,str] %in% c(1, "+"), anno[,start] + zone[1], anno[,end] - zone[2]))
	anno[,end] <- anno[,start] + zone[2] - zone[1]

	##delete duplicates
	if(removeDuplicates)
	{
		anno <- anno[!duplicated(anno),]
	}

	##construct IRanges object
	##FIXME bring in line with earlier checks
	ranges <- IRanges(
		start = anno[,start],
		end = anno[,end]
	)

	##construct GRanges object
	annoZones <- GRanges(
		seqnames = anno[,chr],
		ranges = ranges,
		strand = as(anno[,str], "Rle")
	)

	##get ready to transfer geneIDs
	temp <- DataFrame(anno[,geneID])
	colnames(temp) <- geneID

	##get ready to transfer everything else FIXME test thoroughly
	sel <- !colnames(anno) %in% c(chr,start,end,str,geneID) & !colnames(anno) %in% c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element")
	#1) no redundancies 2) no colnames reserved by GRanges object

	if(any(sel))
	{
		temp2 <- DataFrame(anno[,sel])
		colnames(temp2) <- colnames(anno)[sel]
		temp <- DataFrame(temp, temp2)
	}

	##transfer
	values(annoZones) <- temp

	##output
	annoZones
}
