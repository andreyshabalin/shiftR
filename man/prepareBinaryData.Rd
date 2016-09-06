\name{prepareBinaryData}
\alias{prepareBinaryDataLeft}
\alias{prepareBinaryDataRight}

\title{
	Prepare data for fast circular permutation analysis.
}
\description{
	The algorithm for circular permutation calculation is not symmetric with respect to two datasets and thus required data preprocessing is also different. For simplicity, we call the data sets 'left' and 'right'.
}
\usage{
prepareBinaryDataLeft(set)
prepareBinaryDataRight(set)
}
\arguments{
 \item{set}{
  A 0/1 vector defining selected (genomic) features. The 'left' and 'right' sets must have equal length. The enrichment of their overlap can be assessed with \code{\link{circularPermutationAnalysis}} function.
 }
}

\value{
	Returns objects of class \code{fcpLeft} and \code{fcpRight} respectively.
	The returned objects are used in \code{\link{singleCircularPermutation}} and \code{\link{circularPermutationAnalysis}} functions.
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
 See \code{\link{singleCircularPermutation}} and \code{\link{circularPermutationAnalysis}} function and the respective examples.
}
%\references{}
\examples{
nFeatures = 1e6;

### Generate left and right sets
lset = sample(c(0L,1L), replace = TRUE, nFeatures)
rset = sample(c(0L,1L), replace = TRUE, nFeatures)

# Prepare binary sets:
lbin = prepareBinaryDataLeft(lset)
rbin = prepareBinaryDataRight(rset)

### Check object sizes
# Notice asymetry in binary object sizes

object.size(lset)
object.size(rset)
object.size(lbin)
object.size(rbin)
}
