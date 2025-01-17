% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Asords.R
\name{Asords}
\alias{Asords}
\title{Asymmetric Second Order Rotatable Designs}
\usage{
Asords(v, number_of_pairs, z, type, randomization = FALSE, variance = FALSE)
}
\arguments{
\item{v}{Number of input factors, v(>2)}

\item{number_of_pairs}{Number of pairs of input factors for which asymmetry is required}

\item{z}{A vector of real number and its length equals to number_of_pairs}

\item{type}{Type of central composite design i.e.  ccc or cci. "ccc" is for Central Composite Circumscribed designs and "cci" is for Central Composite Inscribed designs}

\item{randomization}{It is for generating the randomized layout of the design. It takes either TRUE or FALSE  and by default, it is set to FALSE}

\item{variance}{This is for generating the moment matrix and prediction variance of the design based on a second order model. It gives unique prediction variance along with its frequencies. It takes either TRUE or FALSE  and by default, it is set to FALSE}
}
\value{
Asymmetric Second Order Rotatable Designs (ASORDs) for a given v.
}
\description{
This function generates  ASORDs through the orthogonal transformation of central composite designs as per the procedure given by J.S. Mehta and M.N. Das (1968). It would be providing two types of asymmetric designs for a given number of treatments (v). It requires four input parameters viz., v(>2); number_of_pairs(>0); z= vector of real number of length equals to number_of_pairs; type="ccc" or "cci" and randomization=TRUE or FALSE.
}
\examples{

library(AsymmetricSORDs)
Asords(5,2,c(2,3),"ccc",TRUE)

}
\references{
1) J.S. Mehta and M.N. Das (1968)." Asymmetric rotatable designs and orthogonal transformations".

2)M. Hemavathi, Eldho Varghese, Shashi Shekhar & Seema Jaggi (2020)<DOI: 10.1080/02664763.2020.1864817>." Sequential asymmetric third order rotatable designs (SATORDs)".

3) M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI: 10.1080/03610926.2021.1944213>." Theoretical developments in response surface designs: an informative review and further thoughts".
}
