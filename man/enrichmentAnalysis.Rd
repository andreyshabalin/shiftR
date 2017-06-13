\name{enrichmentAnalysis}
\alias{enrichmentAnalysis}
\title{
	Fast enrichment testing via circular permutations on non-binary outcomes.
}
\description{
	This function performs enrichment analysis on
	two sets of matching test statistics.
	The circular permutation scheme accounts for 
	possible local correlation of test statistcs.
	The testing is performed using the quantile thresholds
	provided for each data set.
	
	For every permutation the enrichment is measure with
	Cramer's V coefficient. The maximum/minimum coefficient
	across all considered thresholds is recorded.
	It is then compared with the maximum/minimum coefficient
	observed without permuting the data.
	
	For matching data sets calculated at different genomic locations
	please use \code{\link{matchDatasets}}.
}
\usage{
enrichmentAnalysis(
    pvstats1,
    pvstats2,
    percentiles1 = NULL,
    percentiles2 = NULL,
    npermute,
    margin = 0.05,
    threads = 1)
}
\arguments{
    \item{pvstats1}{
        The vector of statistics for primary data set. \cr
        The statistics must be p-value like, i.e. smaller is better.
    }
    \item{pvstats2}{
        The vector of statistics for auxiliary data set. \cr
        The statistics must be p-value like, i.e. smaller is better.
    }
    \item{percentiles1}{
        These quantile thresholds are used to cut off top results in the primary data set
        for matching with the top results in the auxiliary. \cr
        Can be omitted if the vector \code{pvstats1} is binary.
    }
    \item{percentiles2}{
        Same as \code{percentiles1}, but for the other data set.
    }
    \item{npermute}{
        Number of permutations to perform.
    }
    \item{margin}{
        The minimum offset in the circular permutation to consider. \cr
        Can be a fraction of total number of values or an integer count of values. \cr
        Passed in the call of \code{\link{getOffsetsRandom}} for generation of offsets.
    }
    \item{threads}{
        The number of CPU cores to use for calculations. \cr
        Set to \code{TRUE} to use all cores. \cr
        Multithreading is turned off by default.
    }
}

\value{
	Returns a list with:
	\item{overallPV}{The p-values for the overall test across all thresholds. \cr
	                The p-values are for enrichment, depletion, and two-sided test for both.}
	\item{byThresholdPV}{The p-values for tests for each individual threshold. \cr
	           The p-values provided for enrichment, depletion, and two-sided test.}
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
### Data size
n = 1e5

### Generate vectors of test statistics with local correlation
window = 1000;
pvstats1 = diff(cumsum(runif(n+window)), lag = window)
pvstats2 = diff(cumsum(runif(n+window)), lag = window)

# Add a bit of dependence
pvstats1 = pvstats1 + 0.5 * pvstats2

# test top 0.1, 1, 3, 5, and 10 percent

percentiles1 = c(0.001, 0.01, 0.03, 0.05, 0.1)
percentiles2 = c(0.001, 0.01, 0.03, 0.05, 0.1)

# The offset margin

margin = 0.05;

# Set the number of permutations
# to the maximum

npermute = 1000


enr = enrichmentAnalysis(
        pvstats1,
        pvstats2,
        percentiles1,
        percentiles2,
        npermute,
        margin ,
        threads = 2)

# View the results
enr
}
