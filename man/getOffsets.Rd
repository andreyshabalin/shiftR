\name{getOffsets}
\alias{getOffsetsAll}
\alias{getOffsetsRandom}
\alias{getOffsetsUniform}
\title{
	Generate Random or Uniformly Spaced Permutation Offsets
}
\description{
	This functions generate offsets for permutation analysis with \code{\link{shiftrPermBinary}}. 
	Random, uniformly spaced, and complete sets are available via \code{getOffsetsRandom}, \code{getOffsetsUniform}, and \code{getOffsetsAll} functions respectively.
}
\usage{
getOffsetsRandom(n, npermute, margin = 0.05)
getOffsetsUniform(n, npermute, margin = 0.05)
getOffsetsAll(n, margin)
}
\arguments{
 \item{n}{
  Number of features in the permuted sets.
 }
 \item{npermute}{
  The number of offsets to be generated (number of permutations).
 }
 \item{margin}{
  Offsets by less than \code{margin*n} or more than \code{(1-margin)*n} are not generated.
 }
}

\value{
	Returns a set of permutation offsets for use in \code{\link{shiftrPermBinary}} function. 
	The set of offsets is random for \code{getOffsetsRandom}, uniformly spaced for \code{getOffsetsUniform}, or complete for \code{getOffsetsAll}.

}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
### Number of features, permutations, and margin
n = 100
npermute = 20
margin = 0.1

### Random offsets
getOffsetsRandom(n, npermute, margin = 0.05)

### Uniformly spaced offsets
getOffsetsUniform(n, npermute, margin = 0.05)

### All possible offsets
getOffsetsAll(n, margin)
}
