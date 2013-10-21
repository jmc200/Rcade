\name{exportRcade-methods}
\docType{methods}
\alias{exportRcade}
\alias{exportRcade-methods}
\alias{exportRcade,Rcade-method}

\title{exportRcade and ...}
\description{
Methods for exporting Rcade objects, either to disk or in R.
}

\usage{
exportRcade(x, directory="RcadeOutput", cutoffMode="top", cutoffArg = 1000, justGeneID=FALSE, removeDuplicates="beforeCutoff")
}

\arguments{
  \item{x}{
		An \code{\link[=Rcade-class]{Rcade}} object.
  }
  \item{directory}{
		\code{character} - The directory to export output to.
  }
  \item{cutoffMode}{
		\code{character} - The method to cut off each list (see Details). Must be "all", "top", "B" or "FDR".
  }
  \item{cutoffArg}{
		\code{numeric} - What cutoff to use (see Details).
  }
  \item{justGeneID}{
		\code{logical} - if \code{TRUE}, export only the geneID column. If \code{FALSE}, export all columns.
  }
  \item{removeDuplicates}{
		\code{character} - Should we remove duplicate GeneIDs and, if so, should we do this before or after applying the cutoff? Must be "beforeCutoff", "afterCutoff" or "none". (If removing duplicates then, for each list, the entry with the highest B value is retained.)
  }
}

\section{Usage}{
  \describe{
    \item{}{
      \code{exportRcade(x, directory="RcadeOutput", cutoffMode="top", cutoff = 100, justGeneID=FALSE, removeDuplicates="beforeCutoff", ...)}
    }
  }
}

\details{
	This function exports Rcade output to disk - specifically, it creates the following files:

	\tabular{lll}{
		File:	\tab ChIP:	\tab DE	\cr
		ChIP.csv		\tab	Present (needs log ratio > 0)	\tab	Ignored	\cr
		ChIPonly.csv	\tab	Present (needs log ratio > 0)	\tab	Absent	\cr
		DEandChIP.csv	\tab	Present (needs log ratio > 0)	\tab	Present	\cr
		DownChIP.csv	\tab	Present (needs log ratio > 0)	\tab	Present (logFC < 0)	\cr
		Down.csv		\tab	Ignored	\tab	Present (logFC < 0)	\cr
		DownNoChIP.csv	\tab	Absent	\tab	Present (logFC < 0)	\cr
		Nothing.csv		\tab	Absent	\tab	Absent	\cr
		UpChIP.csv		\tab	Present	(needs log ratio > 0)	\tab	Present (logFC > 0)	\cr
		Up.csv			\tab	Ignored	\tab	Present (logFC > 0)	\cr
		UpNoChIP.csv	\tab	Absent	\tab	Present (logFC > 0)	\cr
	}

	Each file contains genes appropriate to its hypothesis, sorted by descending B value (i.e. ranked from most interesting to least interesting). For example, if you wanted the genes that display DE (either up or down) and also have ChIP signal present, you would look at the top rows of DEandChIP.csv. For genes that have a ChIP signal but explicitly show no DE, use ChIPonly.csv.

	A cutoff is applied to each list, according to the value of \code{cutoffMode}, referring to \code{cutoffArg} if necessary:

	\describe{
		\item{\code{cutoffMode = "all"}}{
			\code{cutoff} ignored, all results written to disk.
		}
		\item{\code{cutoffMode = "top"}}{
			Take the top N genes, where N is specified by \code{cutoffArg}.
		}
		\item{\code{cutoffMode = "B"}}{
			Take all genes with that satisfy \code{B > cutoffArg}, where B is the log-odds.
		}
		\item{\code{cutoffMode = "FDR"}}{
			The expected false positive rate, FPR, and the expected false negative rate, FNR, are calculated using B values.

			The cutoff chosen is the one that maximizes the value of \code{FPR + cutoffArg*FNR}.
		}
	}
}

\examples{
	data(RcadeSTAT1)
	\dontrun{exportRcade(RcadeSTAT1)}
}

\keyword{methods}

