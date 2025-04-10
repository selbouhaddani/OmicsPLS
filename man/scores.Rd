% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OmicsPLS.R
\name{scores}
\alias{scores}
\alias{scores.o2m}
\title{Extract the scores from an O2PLS fit}
\usage{
scores(x, ...)

\method{scores}{o2m}(
  x,
  which_part = c("Xjoint", "Yjoint", "Xorth", "Yorth"),
  subset = 0,
  ...
)
}
\arguments{
\item{x}{Object of class \code{o2m}}

\item{...}{For consistency}

\item{which_part}{character string. One of the following: 'Xjoint', 'Yjoint', 'Xorth' or 'Yorth'.}

\item{subset}{subset of scores vectors to be extracted.}
}
\value{
Scores matrix
}
\description{
This function extracts score matrices from an O2PLS fit
}
\examples{
scores(o2m(scale(-2:2),scale(-2:2*4),1,0,0))

}
\seealso{
\code{\link[OmicsPLS:loadings.o2m]{loadings}}
}
