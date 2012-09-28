\name{constructRcadeTable}
\alias{constructRcadeTable}

\title{
Construct Rcade Table
}
\description{
Most Rcade users will not need to call this function directly. This function constructs a full Rcade table from ChIP and DE data.
}
\usage{
constructRcadeTable(DE, DElookup, chip, annoZone, annoZoneGeneidName, DE.prior=NULL, ChIP.prior=NULL, prior.mode, prior=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DE}{
	\code{data.frame} - DE data (see details section, below)
}
  \item{DElookup}{
	\code{list} - a lookup table specifing the columns of interest in the DE argument.  FIXME - list mandatory columns
}
  \item{chip}{
	\code{data.frame} - ChIP information as ... Columns correspond to samples, and rows should correspond to bins defined by the \code{annoZone} arguments's rows. %%     ~~Describe \code{chip} here~~
}
  \item{annoZone}{
	\code{GRanges} - The genomic bins used in the ChIP-seq analysis. FIXME Metadata must be present.
}
  \item{annoZoneGeneidName}{
	\code{character} - The column in the metadata of \code{annoZone} argument that contains the geneIDs.
}
  \item{DE.prior}{
	FIXME
}
  \item{ChIP.prior}{
	FIXME
}
  \item{prior.mode}{
	FIXME
}
  \item{prior}{
	FIXME
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{RcadeAnalysis}}%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	data(RcadeSTAT1)

	dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

	DE <- getDE(RcadeSTAT1)
	DElookup <- list(GeneID="ENSG", logFC="logFC", B="B",
		"Genes.Location", "Symbol")

	chip <- getChIP(RcadeSTAT1)
	annoZone <- getChIP(RcadeSTAT1, what="annoZones")

	x <- constructRcadeTable(DE, DElookup, chip, annoZone, annoZoneGeneidName="ENSG", prior.mode="assumeIndependent")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
