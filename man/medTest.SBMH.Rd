\name{medTest.SBMH}
\alias{medTest.SBMH}
\title{
Test for multiple mediators which uses the Bogomolov and Heller approach
}
\description{
Implements Heller/Sampson method of testing multiple biological mediators
simultaneously, controlling either FWER or FDR.
}
\usage{
medTest.SBMH(pEM,pMY,MCP.type="FWER",t1=0.05,t2=0.05,lambda=0)
}
\arguments{
  \item{pEM}{
    Vector of size m (where m = number of mediators). Entries are the p-values for the E,M_j relationship.
  }
  \item{pMY}{
  Vector of size m (where m = number of mediators). Entries are the p-values for the M_j,Y|E relationship.
  }
  \item{MCP.type}{
  Multiple comparison procedure - either "FWER" or "FDR".
  }
  \item{t1}{
  Threshold for determining the cutoff to be one of the top S_1 E/M_j relationships.
  }
  \item{t2}{
  Threshold for determining the cutoff to be one of the top S_2 M_j/Y relationships.
  }
  \item{lambda}{
  Threshold for estimating the proportion of false positives.
  }
}
\details{
  See Heller/Sampson paper for a more in-depth description of this method.
}
\value{
  m x 1 matrix of either p-values (if MCP.type = "FWER") or q-values (if MCP.type = "FDR").
}
\references{
}
\author{
  Ruth Heller, Joshua N. Sampson
}

\examples{

}
\keyword{htest}
