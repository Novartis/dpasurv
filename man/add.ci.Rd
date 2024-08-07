% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{add.ci}
\alias{add.ci}
\title{Calculate bootstrap confidence bands for effect of interest}
\usage{
add.ci(object, alpha)
}
\arguments{
\item{object}{object of class "effect"}

\item{alpha}{the confidence level}
}
\value{
object of class "effect" with updated confidence interval corresponding to alpha
}
\description{
Calculate bootstrap confidence bands for effect of interest
}
\examples{
library(dpasurv)

data(simdata)

# Perform dynamic path analysis:
# We set boot.n=30 for the example to run fast, should be set large enough
# so that results don't change meaningfully for different seeds.
s <- dpa(Surv(start,stop,event)~M+x, list(M~x), id="subject", data=simdata, boot.n=30)

# Calculate cumulative direct effect (which calculates a CI with alpha = 0.05 by default):
direct <- effect(x ~ outcome, s)

# update confidence interval for a new alpha (this overwrites the 0.05 CI already calculated above)
direct <- add.ci(direct, alpha=0.10)
}
