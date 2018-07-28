\name{NavyAdenoma}
\docType{data}
\alias{NavyAdenoma}

\title{Data from the Sinha et al. 1999 and Boca et al. 2014 studies
  concerning the Navy Colorectal Adenoma case-control study.}
\description{The data consists of a dietary exposures of interest (daily
  intake of fish in grams/day), the outcome of interest (colorectal
  adenoma status), and 149 serum metabolites which are possible
  mediators between either of the exposures and the outcome. Four other
  additional covariates of interest are included (body mass index,
  gender, age, current smoking status.) }

\usage{data(NavyAdenoma)}
\format{
  Data frame with 258 rows, representing 129 cases and 129 controls in
  the Navy Colorectal Adenoma study. There are 155 columns, ordered as
  follows: fish intake (\emph{Fish}), body mass index (\emph{BMI}),
  gender (\emph{Female} - coded as 0 for male, 1 for female), age in
  years (\emph{Age}), current smoking status (\emph{Smoking} - coded as
  0 for non-smoker, 1 for current smoker), 149 serum metabolites (values
  were previously batch normalized and log transformed), and
  case-control 
  status for adenoma (\emph{Adenoma} - coded as 0 for control, 1 for case.)
}

\references{
  Sinha, R, Chow WH, Kulldorff M, Denobile J, Butler J, Garcia-Closas M,
  Weil R, Hoover RN, Rothman N. Well-done, grilled red meat
  increases the risk of colorectal adenomas. \emph{Cancer Res.}, 1999, 59(17),
  4320-4324.
  
  Boca SM, Sinha R, Cross AJ, Moore SC, Sampson JN. 
  Testing multiple biological mediators simultaneously.
  \emph{Bioinformatics}, 2014, 30(2), 214-220.
}
\seealso{
  \code{medTest}
}

\keyword{datasets}
