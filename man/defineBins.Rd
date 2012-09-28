\name{defineBins}
\alias{defineBins}

\title{Define Bins}

\description{
Defines bins about the 5' end of certain features of interest - these features are usually transcripts.
}

\usage{
defineBins(anno, zone, geneID="ensembl_gene_id", removeDuplicates=TRUE)
}

\arguments{
  \item{anno}{
	\code{data.frame} (or, an object that can be coerced to a data.frame, such as a GRanges) - Annotation information, corresponding to features of interest (usually transcripts). Only the 5' end of each object is used.

	Rcade expects the following column names: \code{chr, start, end, str}. These correspond to chromosome name, start co-ordinate, end co-ordinate and strand. Additionally, there must be another column specifying a gene ID, specified by the \code{geneID} argument.
}
  \item{zone}{\code{integer} - must be a length 2 vector of form \code{c(relative.start, relative.end)}. For example, \code{zone = c(-10,100)} will produce bins that start 10bp 5' of each transcript's TSS and end 100bp 3' of it.}

  \item{geneID}{\code{character} or \code{integer} - The column in \code{anno} that contains a geneID (or some other feature ID).}

  \item{removeDuplicates}{\code{logical} - If \code{TRUE}, then any rows that share the same geneID and genomic location as another row will be removed (even if any of the other columns are different).}
}

\details{
The \code{defineBins} function is useful when ChIP-seq bins are defined about ... . In particular, biomaRt data can be fed into this function directly. FIXME See vignette.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A \code{\link{GRanges}} object, corresponding to genomic bins. This output can be used as the \code{ChIPannoZones} argument in \code{\link{RcadeAnalysis}}.
}

\author{Jonathan Cairns}

\seealso{
\code{\link{RcadeAnalysis}}
}
\examples{
\dontshow{
	#read example annotation - truncated, do NOT use this in your analysis!
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")
	anno <- read.csv(file.path(dir, "anno.csv"))
}

\dontrun{	##acquire annotation from biomaRt
	library(biomaRt)
	anno <- getBM(
		attributes= c("ensembl_gene_id", "chromosome_name",
			"transcript_start", "transcript_end", "strand"),
		mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	)
}
#define bins about the annotation
anno <- anno[order(anno$chromosome_name),]
colnames(anno) <- c("ENSG","chr","start","end","str")
ChIPannoZones <- defineBins(anno, c(-1500, 1500), geneID = "ENSG")

}

