\name{diffCountsBaySeq}
\alias{diffCountsBaySeq}

\title{
Differential Counts wrapper - BaySeq
}
\description{
Most Rcade users will not need to call this function directly. A function that provides a wrapper for the methods in the BaySeq package. 
}
\usage{
diffCountsBaySeq(counts, targets, annoZones, cl = NULL, getLibsizesArgs = list(estimationType = "quantile", quantile = 0.75), getPriors.NBArgs = list(), getLikelihoods.NBArgs = list(), libsizes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
		Counts from \code{\link{countReads}}%%     ~~Describe \code{counts} here~~
	}
  \item{targets}{
		\code{Data.frame} - Information about the ChIP data files. Mandatory column names are: "fileid", "sampleid", "factor", "filepath".
	}
  \item{annoZones}{
		\code{GRanges} specifying the bins of interest, with a column in the metadata for the geneID.%%FIXME copy from RcadeAnalysis.Rd
	}
  \item{cl}{
	\code{cluster} from \code{\link[parallel]{makeCluster}} in the parallel package.
	}
  \item{getLibsizesArgs}{
	\code{List} - Arguments to be passed to the getLibsizes function. If a \code{libsizes} column is present in the targets file, then these arguments are ignored. %FIXME...

	\code{getLibsizesArgs$cD} is always ignored.

	See \code{\link[baySeq]{getLibsizes}} for a list of arguments. %FIXME link across packages
	}
  \item{getPriors.NBArgs}{
	See \code{\link[baySeq]{getPriors}} for a list of arguments. %FIXME link across packages

	\code{getPriors.NBArgs$cD} and \code{getPriors.NBArgs$cl} are always ignored.
	}
  \item{getLikelihoods.NBArgs}{
	See \code{\link[baySeq]{getLikelihoods}} for a list of arguments. %FIXME link across packages

	\code{getLikelihoods.NBArgs$cD} and \code{getLikelihoods.NBArgs$cl} are always ignored.
	}
  \item{libsizes}{
		Library sizes FIXME
	}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	\code{data.frame} containing differential count information. %FIXME
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Jonathan Cairns}
\note{
%%  ~~further notes~~
}

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

	counts <- countReads(ChIPannoZones, targets, fileDir = dir)

	x <- diffCountsBaySeq(counts, targets, ChIPannoZones)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
