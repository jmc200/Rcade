\name{Rcade-class}
\docType{class}
\alias{Rcade-class}

\alias{plotBBB}
\alias{plotMM}
\alias{plotPCA}

\alias{plotBBB,Rcade-method} %FIXME document all of these
\alias{plotMM,Rcade-method}
\alias{plotPCA,Rcade-method}

\alias{getDE}
\alias{getChIP}
\alias{getRcade}

\alias{getDE,Rcade-method}
\alias{getChIP,Rcade-method}
\alias{getRcade,Rcade-method}

\alias{getDE<-}
\alias{getChIP<-}
\alias{getRcade<-}

\alias{getDE<-,Rcade-method}
\alias{getChIP<-,Rcade-method}
\alias{getRcade<-,Rcade-method}

\title{Rcade Class}
\description{
The main class in Rcade. This class contains data pertaining to any relevant DE experiments, ChIP-seq experiments, and Rcade output from linking the previous two.

Objects of this class are typically created with the \code{\link{RcadeAnalysis}} function.
}

\section{Plotting methods}{
  \describe{
    \item{}{
      \code{plotPCA(x, ...)}: Perform PCA analysis on the ChIP-seq data and plot the results. 
	}
    \item{}{
      \code{plotMM(x, ...)}: Plot ChIP log-ratios against DE log-ratios. \code{...} arguments are passed to \code{\link{plot}}. FIXME 
    }
    \item{}{
      \code{plotBBB(x, ...)}: (NB: Requires the CRAN package \pkg{rgl}.) 3D plot comparing log-odds values for ChIP, DE and combined ChIP & DE. \code{...} arguments are passed to \code{\link{plot}}.
    }
  }
}

\section{Accessors}{
  \describe{
    \item{}{
      \code{getDE(x, what="summary")}: Get DE information. \code{what} can be: \code{"summary"} for the DE analysis, \code{"prior"} for the prior probability/probabilities of DE presence.
    }
    \item{}{
      \code{getChIP(x, what="summary")}: Get ChIP analysis information. \code{what} can be: \code{"summary"} for the analysis, \code{"counts"} for the raw counts, \code{"annoZones"} for the bins used in the analysis, \code{"prior"} for the prior probability/probabilities of ChIP signal presence, or \code{"targets"} for the targets file.
    }
    \item{}{
      \code{getRcade(x)}: Get the Rcade table - i.e. combined DE/ChIP information.
    }
  }
}

\author{Jonathan Cairns}

\seealso{\code{\link{RcadeAnalysis}}}

\examples{
	data(RcadeSTAT1)
	RcadeSTAT1

	x <- getChIP(RcadeSTAT1)
	y <- getDE(RcadeSTAT1)
	z <- getRcade(RcadeSTAT1)

	plotMM(RcadeSTAT1)
	plotPCA(RcadeSTAT1)
	library(rgl) ##required for plotBBB
	plotBBB(RcadeSTAT1)
}

\references{
	NA
}
\keyword{misc}
