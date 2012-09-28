\name{countReads}
\alias{countReads}

\title{
Count Reads
}
\description{
Most Rcade users will not need to call this function directly. Given targets information linking to bam files, count the reads that lie in defined bins.
}
\usage{
countReads(annoZone, targets, fileDir=NULL, dontCheckTargets=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annoZone}{
		\code{GRanges} - The bins to be used when counting reads.
}
  \item{targets}{
		\code{data.frame} - Targets file (see vignette)
}
  \item{fileDir}{
		\code{character} - The directory in which the raw ChIP-seq data files are kept.
}
  \item{dontCheckTargets}{
		\code{logical} - If \code{TRUE}, the targets file is not checked for consistency/appropriate field names. This should not be changed for Rcade purposes, but may be useful if you wish to obtain bin counts for some other purpose. Make sure relevant column names are lower case. Use at your own risk!
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Matrix of read counts, with columns corresponding to samples and rows corresponding to bins.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Jonathan Cairns}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{RcadeAnalysis}}
}
\examples{
	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	targets <- read.csv(file.path(dir, "targets.csv"), as.is = TRUE)

	anno <- read.csv(file.path(dir, "anno.csv"))
	anno <- anno[order(anno$chromosome_name),]
	colnames(anno) <- c("ENSG","chr","start","end","str")
	ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")

	x <- countReads(ChIPannoZones, targets, fileDir = dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
