\name{medTest}
\alias{medTest}
\title{
Permutation test for multiple mediators which uses a joint correction 
}
\description{
Implements Boca et al. 2014 method of testing multiple biological mediators
simultaneously, which uses a permutation test with a joint correction based
on the maximal test statistics.
}
\usage{
medTest(E, M, Y, Z = NULL, useWeightsZ = TRUE, nperm = 100, w = 1)
}
\arguments{
  \item{E}{
    Exposure: numerical vector of length n or numerical n x 1 matrix,
    where n is the number of samples or subjects.
  }
  \item{M}{
    Mediators: numerical matrix of dimension n x p, where p is the
    number of mediators.
  }
  \item{Y}{
    Outcome: numerical vector of length n or numerical n x 1 matrix
  }
  \item{Z}{
    Additional covariates: NULL or numerical matrix having n rows. If it
    is not NULL, then the exposure, mediators, and outcome will all be
    initially regressed on Z, with the residuals being used in the
    mediation analysis.
  }
  \item{nperm}{
    Number of permutations performed. One of the permutations always
    consists of the observed data.
  }
  \item{w}{
    Weight assigned to each subject or sample for the E-M associations.
    Useful for case-control sampling.
  }
  \item{useWeightsZ}{
    If TRUE, use the weights in w for the initial regression of E and M
    on Z, in the case where Z is not null. If FALSE, use equal weights.
  }
}
\details{
  See Boca et al. 2014 paper for a more in-depth description of this method.
}
\value{
  Matrix with p rows and 2 columns, where the first column (labeled "S")
  gives the S-statistics and the second column (labeled "p")
  gives the p-values for each mediator, where the mediators correspond to the 
  rows. For each mediator, the S-statistic is the absolute value of the product
  of the correlation between the exposure and the mediator and the
  conditional (or partial) correlation between the mediator and the
  outcome conditional on the exposure.  Thus, the larger the S-statistic
  is, the smaller the p-value is, and the more likely it is
  that the mediator is statistically significant.
}
\references{
  Boca SM, Sinha R, Cross AJ, Moore SC, Sampson JN. (2014) 
  Testing multiple biological mediators simultaneously.
  \emph{Bioinformatics}, 30(2), 214-220.
}
\author{
  Simina M. Boca, Joshua N. Sampson
}

\seealso{
  \code{NavyAdenoma}
}     
\examples{
##More details on this example are given in the accompanying vignette.

data(NavyAdenoma)

##the exposure of interest is the daily intake of fish
##the possible mediators are 149 normalized serum metabolite values
##the outcome of interest is colorectal adenoma case-control status

##get weights needed due to retrospective sampling, to use for all regressions
##that do not have adenoma status as the outcome
prev <- 0.228
p <- sum(NavyAdenoma$Adenoma==1)/nrow(NavyAdenoma)
w <- rep(NA, nrow(NavyAdenoma))
w[NavyAdenoma$Adenoma == 1] <- prev/p
w[NavyAdenoma$Adenoma == 0] <- (1-prev)/(1-p)

##perform the test of mediation, adjusting for the covariates
##BMI, gender, age, and current smoking status
##(only running 10 permutations in this example - would need to
##perform more to get an accurate p-value)
set.seed(840218)
medsFish <- medTest(E=NavyAdenoma$Fish,
                    M=NavyAdenoma[, 6:154],
                    Y=NavyAdenoma$Adenoma,
                    Z=NavyAdenoma[, 2:5],
                    nperm=10, w=w)

##get metabolite with lowest p-value
medsFish[which.min(medsFish[,"p"]),]
colnames(NavyAdenoma[,-(1:5)])[which.min(medsFish[,"p"])]
}
\keyword{htest}
