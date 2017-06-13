\name{cramerV}
\alias{cramerV}
\title{
	Calculate Cramer's V (phi) coefficient
}
\description{
	This functions calculates Cramer's V coefficient
	for overlap of two binary data sets.
}
\usage{
cramerV(sum12, sum1, sum2, len)
}
\arguments{
    \item{len}{
        Total number of elements in each data set.
    }
    \item{sum1}{
        Number of active features in data set 1.
    }
    \item{sum2}{
        Number of active features in data set 1.
    }
    \item{sum12}{
        Number of jointly active features in both data sets.
    }    
}

\value{
	Returns the Cramer's V coefficient.
}
\note{
    The parameters can be single values or vectors.
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
# Zero score for perfect independence
cramerV(100,10000,10000,1000000)

# Positive score for increased overlap
cramerV(150,10000,10000,1000000)

# Negative score for decreased overlap
cramerV( 50,10000,10000,1000000)

# We can input a vector for sum12
cramerV(99:101,10000,10000,1000000)
}
