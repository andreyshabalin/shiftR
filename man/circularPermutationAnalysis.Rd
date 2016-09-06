\name{circularPermutationAnalysis}
\alias{circularPermutationAnalysis}
\title{
	Fast Enrichment Testing via Circular Permutations
}
\description{
	This function performs very fast feature enrichment permutation testing between two binary data sets.
	Circular permutations are used instead of simple permutations to preserve local dependence of test statistics.
	The input data sets must be prepared with \code{\link{prepareBinaryDataLeft}} and \code{\link{prepareBinaryDataRight}} functions.
}
\usage{
circularPermutationAnalysis(left,
                          right,
                          npermute,
                          margin = 0.05,
                          alsoDoFisher = TRUE,
                          returnPermOverlaps = FALSE)
}
\arguments{
 \item{left}{
  Feature set prepared with \code{\link{prepareBinaryDataLeft}} function.
 }
 \item{right}{
  Feature set prepared with \code{\link{prepareBinaryDataRight}} function.
 }
 \item{npermute}{
  Number of permutations to perform.
 }
 \item{margin}{
  Circular permutation by less than \code{margin} of the whole feature set are not considered.
 }
 \item{alsoDoFisher}{
  If \code{TRUE}, also perform Fisher exact test (via \code{\link[stats]{fisher.test}}).
 }
 \item{returnPermOverlaps}{
  If \code{TRUE} return overlap counts under all tested permutations. 
 }
}

\value{
	Returns a list with:
	\item{nfeatures}{Number of features in input data sets.}
	\item{lfeatures}{Number of active features in the left data set.}
	\item{rfeatures}{Number of active features in the right data set.}
	\item{overlap}{Number of features simultaneously active in both data sets.}
	\item{overlapUnderNull}{Expected value of \code{overlap} if input data sets were independent.}
	\item{enrichment}{ Enrichment ratio, equal to \code{overlap / overlapUnderNull}}
	\item{permPVenrich}{ Permutation p-value for enrichment (one-sided).}
	\item{permPVdeplete}{ Permutation p-value for depletion (one-sided).}
	\item{permPV}{ Permutation p-value for depletion (two-sided).}
	\item{permZ}{Permutation z-statistic, calculated by fitting normal distribution to the \code{overlap} values under permutations. \cr
					 Positive values indicate enrichment.}
	\item{fisherTest}{Fisher exact test, as output by \code{\link[stats]{fisher.test}}}
	\item{fisherMat}{Input 2x2 matrix for Fisher exact test.}
	\item{overlapsPerm}{Vector of length \code{npermute} with \code{overlap} values under permutations.}
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
 This function essentially involves \code{npermute} calls of \code{\link{singleCircularPermutation}} function.
}
\examples{
### Number of features
nf = 1e6;

### Generate data sets
# The vector of a few common active feature to create dependence
common = sample(c(0L,1L), size = nFeatures, replace = TRUE, prob = c(0.999,0.001))
# Left and right data sets with the common active features
lset = sample(c(0L,1L), size = nf, replace = TRUE, prob = c(0.8,0.2)) | common
rset = sample(c(0L,1L), size = nf, replace = TRUE, prob = c(0.8,0.2)) | common

# Prepare binary sets:
left = prepareBinaryDataLeft(lset)
right = prepareBinaryDataRight(rset)


z = circularPermutationAnalysis(left, right, npermute = 10000)

show(z)
}
