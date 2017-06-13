\name{shiftrPrepare}
\alias{shiftrPrepareLeft}
\alias{shiftrPrepareRight}

\title{
	Prepare data for fast circular permutation analysis.
}
\description{
    The concept of circular permutations
    is symmetric with respect to the input data sets.
	The algorithm for circular permutation calculation is, however,
	not symmetric with respect to two datasets and 
	thus the required data preprocessing is also different.
	For simplicity, we call the data sets 'left' and 'right'.
}
\usage{
shiftrPrepareLeft(set)
shiftrPrepareRight(set)
}
\arguments{
    \item{set}{
        A 0/1 vector defining selected (genomic) features. 
        The 'left' and 'right' sets must have equal length. 
        The enrichment of their overlap can be assessed w
        ith \code{\link{shiftrPermBinary}} function.
    }
}

\value{
	Returns objects of class \code{fcpLeft} and \code{fcpRight} respectively.
	The returned objects are used in \code{\link{singlePermutation}} and \code{\link{shiftrPermBinary}} functions.
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
 See code{\link{shiftrPermBinary}} function and the respective example.
}
%\references{}
\examples{
### Number of features
nf = 1e6

### Generate left and right sets
lset = sample(c(0L,1L), size = nf, replace = TRUE)
rset = sample(c(0L,1L), size = nf, replace = TRUE)

# Prepare binary sets:
lbin = shiftrPrepareLeft(lset)
rbin = shiftrPrepareRight(rset)

### Check object sizes
# Notice asymetry in binary object sizes

object.size(lset)
object.size(rset)
object.size(lbin)
object.size(rbin)
}
