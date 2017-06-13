\name{matchDatasets}
\alias{matchDatasets}
\title{
	Match two data sets by location
}
\description{
    The goal of this function is to match records in the data sets
    for subsequent enrichment analysis.

	For each record in the primary data set (\code{data1})
	it finds the record in the auxiliary data set (\code{data1})
	which overlap with it or lie within the flanking distance (\code{flank}).
	If multiple such auxiliary record are found, 
	we select the one with the center closest to
	the center of the primary record.
	If no such record is available, no matching is made for the primary record.
}
\usage{
matchDatasets(  data1, data2, flank = 0)
}
\arguments{
    \item{data1}{
        A data frame with the primary data set, must have at least 4 columns:
        \enumerate{
            \item Chromosome name.
            \item Start position.
            \item End position.
            \item P-value or test statistic.
            \item Optional additional columns.
        }
    }
    \item{data2}{
        A data frame with the auxiliary data set. \cr
        Must satisfy the same format criteria as the primary data set.
    }
    \item{flank}{
        Allowed distance between matched records. \cr
        Set to zero to require overlap.
    }
}

\value{
	Returns a list with three elements.
	\item{data1}{The primary data sets with unmatched records excluded.}
	\item{data2}{The auxiliary data set records matching those in \code{data1} above.}
}
\note{
    For a technical reason, the chromosome positions are assumed to be no greater than \code{1e9}.
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{

data1 = read.csv(text =
"chr,start,end,stat
chr1,100,200,1
chr1,150,250,2
chr1,200,300,3
chr1,300,400,4
chr1,997,997,5
chr1,998,998,6
chr1,999,999,7")

data2 = read.csv(text =
"chr,start,end,stat
chr1,130,130,1
chr1,140,140,2
chr1,165,165,3
chr1,200,200,4
chr1,240,240,5
chr1,340,340,6
chr1,350,350,7
chr1,360,360,8
chr1,900,900,9")

# Match data sets exactly
matchDatasets(data1, data2, 0)

# Match data sets with a flank
# last records are now matched
matchDatasets(data1, data2, 100)
}
