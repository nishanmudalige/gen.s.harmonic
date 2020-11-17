\name{y.q.x}
\alias{y.q.x}
\title{A function to calculate the generalized spherical harmonics}
\description{A function to calculate the generalized spherical harmonics as defined by Dai and Xu in "Approximation Theory and Harmonic Analysis on Spheres and Balls".}
\usage{y.q.x(x, q, normalize.by.volume)}
\arguments{
  \item{x}{A vector of the form x = (x_1, ..., x_p) where p>=3 and each x_j is a real number.}
  \item{q}{A vector of the form q = (q_1, ..., q_p) where p>=3 and each q_j is a non-negative integer.}
  \item{normalize.by.volume}{Logical argument. Whether to normalize by the volume of the hyper-sphere. If TRUE, then the value returned with be normalized by the volume element.}
}
\value{The generalized spherical harmonic evaluated for x evaluated for q.}
\references{Dai, F. and Xu, Y. (2013), Approximation Theory and Harmonic Analysis on Spheres and Balls, Springer Science & Business Media.}
\author{Nishan Mudalige}
\note{This research was supported by a grant from the Google Cloud Platform.}
\examples{
## Example 1
# q = c(1,0,0)
# x = c(1,1,1)

# y.q.x(x, q, F)
# y.q.x(x, q, T)


## Example 2
# q = c(1,0,0)
# x = c(1/sqrt(3),1/sqrt(3),1/sqrt(3))

# y.q.x(x, q, F)
# y.q.x(x, q, T)


## Example 3
# q = c(1,0,0,0)
# x = c(1/sqrt(8),2/sqrt(8),3/sqrt(8),4/sqrt(8))

# y.q.x(x, q, F)
# y.q.x(x, q, T)



## Example 4
# q = c(0,1,0,0)
# x = c(1/sqrt(8),2/sqrt(8),3/sqrt(8),4/sqrt(8))

# y.q.x(x, q, F)
# y.q.x(x, q, T)
}
