% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pred.var.R
\name{Pred.var}
\alias{Pred.var}
\title{Function for generating the moment matrix and variance of the predicted response}
\usage{
Pred.var(matrix)
}
\arguments{
\item{matrix}{Design matrix with the coefficients of the corresponding input factors}
}
\value{
The moment matrix and the prediction variance for a given design based on a second-order model
It gives unique prediction variance along with its frequencies.
}
\description{
This function generates the moment matrix and variance of the predicted response for a given design based on a second-order model, for measuring the rotatability of the design. The input should be the specified form of a design matrix with the coefficients of the corresponding input factors. A minimum number of centre points is to be used to ensure the non-singularity of X`X.
}
\examples{
 \dontrun{
library(AsymmetricSORDs)
Pred.var(matrix)
}
}
\references{
1) G.E.P. Box and K.B. Wilson (1951).' On the experimental attainment of optimum conditions'.

2) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>.' Theoretical developments in response surface designs: an informative review and further thoughts'.
}
