% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCD_coded.R
\name{CCD_coded}
\alias{CCD_coded}
\title{Central Composite Designs (CCD) with coded levels}
\usage{
CCD_coded(v, type, randomization = FALSE, variance = FALSE)
}
\arguments{
\item{v}{Number of input factors, v(>2)}

\item{type}{Type of central composite design i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs}

\item{randomization}{It is for generating the randomized layout of the design. It takes either TRUE or FALSE  and by default, it is set to FALSE}

\item{variance}{This is for generating the moment matrix and prediction variance of the design based on a second order model. It gives unique prediction variance along with its frequencies. It takes either TRUE or FALSE  and by default, it is set to FALSE}
}
\value{
Central Composite Designs (CCD) for a given number of input factors (v) with coded levels
}
\description{
This function generates Central Composite Designs (CCD) with coded levels for a given number of input factors (v). The CCD constitute  combinations of factorial points, axial points and center points. Three types of CCD can be generated using this function  i.e.  ccc or cci or ccf. "ccc" is for Central Composite Circumscribed designs, "cci" is for Central Composite Inscribed designs and "ccf" is for Central Composite Face Centered designs. It gives the randomized layout of the design along with the moment  matrix and prediction variance.
}
\note{
Here, the factorial portion consists of  2^v (full factorial) combinations and there is no upper limit for the number of input factors, v (>2). To get a CCD with smaller runs,  one may use fractional factorial (of resolution V) in place of full factorial.
}
\examples{

library(AsymmetricSORDs)
CCD_coded(5,'ccc',FALSE,FALSE)
CCD_coded(6,"cci",FALSE,FALSE)

}
\references{
1) G.E.P. Box and K.B. Wilson (1951)." On the experimental attainment of optimum conditions".

2) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>. "Theoretical developments in response surface designs: an informative review and further thoughts".
}
