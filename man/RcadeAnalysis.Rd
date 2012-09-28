\name{RcadeAnalysis}
\alias{RcadeAnalysis}

\title{Rcade Analysis}
\description{
The main function in Rcade - reads in DE information, processes ChIP data from raw .bam files, and then combines the two to form an Rcade object.
}
\usage{
RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName, ChIPtargets, ChIPfileDir, cl, DElookup, DE.prior = NULL, prior.mode="assumeIndependent", prior=NULL, ...)
}
\arguments{
  \item{DE}{\code{Data.frame} - DE summary information for genes of interest. For example, output from limma. EITHER \code{DE} must have column names "geneID", "logfc" and "B" (case insensitive) OR you should specify DElookup.}
  \item{ChIPannoZones}{\code{GRanges} specifying the bins of interest, with a column in the metadata for the geneID.} %FIXME function
  \item{annoZoneGeneidName}{\code{character} - the name of the column in \code{ChIPannoZones}'s metadata corresponding to geneID.}
  \item{ChIPtargets}{\code{Data.frame} - Information about the ChIP data files. Mandatory column names are: "fileid", "sampleid", "factor", "filepath".}
  \item{ChIPfileDir}{\code{character} - Directory, within which "filepath" of \code{ChIPtargets} is evaluated.}
  %\item{ChIPshift}{\code{Integer} - Reads on the forward strand are moved to the right by this number of base pairs, and reads on the negative strand are moved to the left}
  \item{cl}{A \code{cluster} from \code{\link[parallel]{makeCluster}} in the parallel package.}
  \item{DElookup}{\code{list} - lookup table of form \code{list(RcadeField1 = DEcolumn1, RcadeField2 = DEcolumn2, ...)}. If you don't specify this argument, then Rcade will try to find the mandatory fields automatically but will not keep any of the other information in its output.}
  \item{DE.prior}{\code{numeric} - Currently ignored.}%FIXME The prior probability of DE for each GeneID. Either a scalar, or a vector where the Nth element corresponds to the Nth row of the DE argument. }
  \item{prior.mode}{The method used to create prior probabilities in the Rcade table. "assumeIndependent" is currently the only option available.} %FIXME
  \item{prior}{Currently ignored.} %The prior to be used.} %FIXME
  \item{...}{Additional arguments.}
}
\details{
	This is the main analysis function in Rcade. The user should specify information relating to the DE and ChIP data for the experiment in question. Rcade will process these data and rank genes by the combined DE and ChIP strength.
}

\value{
	An Rcade object.
}

\author{Jonathan Cairns}

\seealso{\code{\link{RcadeAnalysis}}}

\examples{
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	DE <- read.csv(file.path(dir, "DE.csv"))
	DElookup <- list(GeneID="ENSG", logFC="logFC", B="B",
"Genes.Location", "Symbol")

	targets <- read.csv(file.path(dir, "targets.csv"), as.is = TRUE)
	
	anno <- read.csv(file.path(dir, "anno.csv"))
	anno <- anno[order(anno$chromosome_name),]
	colnames(anno) <- c("ENSG","chr","start","end","str")
	ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")

	Rcade <- RcadeAnalysis(DE, ChIPannoZones, ChIPtargets=targets, ChIPfileDir = dir,
	DElookup=DElookup)
}

%\references{
%	NA
%}
\keyword{misc}
