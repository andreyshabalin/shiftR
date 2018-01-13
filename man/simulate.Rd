\name{simulate}
\alias{simulateNumeric}
\alias{simulateBinary}
\alias{simulatePValues}
\title{
    Generate Artificial Data for Tests and Illustrations
}
\description{
    These functions generate two artificial data sets
    with local dependence of observations.
}
\usage{
simulateNumeric(n, corWithin, corAcross = 0)
simulateBinary(n, corWithin, corAcross = 0)
}
\arguments{
    \item{n}{
        Total number of elements in each data set.
    }
    \item{corWithin}{
        Correlation of adjacent observations within each data set.
    }
    \item{corAcross}{
        Correlation of observations across data sets.
    }
}

\value{
    Returns the Cramer's V coefficient.
}
\note{
    The \code{simulateNumeric} function generates two data sets with elements
    having standard normal distribution. 
    
    The \code{simulateBinary} function generates data sets with 0/1 values
    by thresholding the numeric data sets from \code{simulateNumeric}.
    
    The \code{simulatePValues} function generates data sets of p-values
    by applying \code{\link[stats:Normal]{pnorm}} to the data sets 
    from \code{simulateNumeric}.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\examples{
n = 100000
sim = simulateNumeric(n, 0.5, 0.3)

# Means should be close to 0 (zero)
mean(sim$data1)
mean(sim$data2)

# Variances should be close to 1
var(sim$data1)
var(sim$data2)

# Correlation of adjacent observations
# should be close to 0.5
cor(sim$data1[-1], sim$data1[-n])
cor(sim$data2[-1], sim$data2[-n])

# Correlation between data sets 
# should be close to 0.3
cor(sim$data1, sim$data2)
}
