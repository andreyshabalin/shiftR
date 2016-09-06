\name{singleCircularPermutation}
\alias{singleCircularPermutation}
\title{
	Count feature overlap under permutation (offset).
}
\description{
	This function performs fast feature overlap count under a circular permutation.
	The input data sets must be prepared with \code{\link{prepareBinaryDataLeft}} and \code{\link{prepareBinaryDataRight}} functions.
}
\usage{
singleCircularPermutation(left, right, offset)
}
\arguments{
 \item{left}{
  Feature set prepared with \code{\link{prepareBinaryDataLeft}} function.
 }
 \item{right}{
  Feature set prepared with \code{\link{prepareBinaryDataRight}} function.
 }
 \item{offset}{
  Offset of one feature set relative to another. See the example below for clarity. \cr
  Zero indicate no offset, i.e. simply count feature overlap.
 }
}

\value{
	Returns count of feature overlap under a circular permutation (offset).
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
### Number of features
nf = 1e6

### Generate left and right sets
lset = sample(c(0L,1L), size = nf, replace = TRUE)
rset = sample(c(0L,1L), size = nf, replace = TRUE) | lset

# Prepare binary sets:
lbin = prepareBinaryDataLeft(lset)
rbin = prepareBinaryDataRight(rset)

### count feature overlap
# R calculations
overlapS = sum(lset & rset)
# Binary calculations
overlapF = singleCircularPermutation(lbin, rbin, 0)

message("Feature overlap: ",
        overlapS, " / ", overlapF,
        " (slow/fast count)")
stopifnot( overlapS == overlapF )


### Count overlap with offset
offset = 2017
# R calculations
overlapOS = sum(lset[ c((offset+1):nf, 1:offset)] & rset)
# Binary calculations
overlapOF = singleCircularPermutation(lbin, rbin, offset)

message("Feature overlap at offset: ",
        overlapOS, " / ", overlapOF,
        " (slow/fast count)")
stopifnot( overlapOS == overlapOF )
}
