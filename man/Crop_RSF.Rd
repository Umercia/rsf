% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_SF.R
\name{Crop_RSF}
\alias{Crop_RSF}
\title{Crop_RSF Function}
\usage{
Crop_RSF(rsf, Xmax, Xmin, Ymax, Ymin)
}
\arguments{
\item{rsf:}{rsf object already loaded into memory}

\item{Xmax:}{X maximum limit}

\item{Xmin:}{X minimum limit}

\item{Ymax:}{Y maximum limit}

\item{Ymin:}{Y minimum limit}
}
\description{
Crop an *.rsf file around a rectangle define by Xmax, Xmin, Ymax, Ymin.
}
\examples{
load into memory an rsf file into the object rsf_166:
rsf_166 <- Read_RSF("RSF-windresource-CFD_Aldermyrberget 166m.rsf")
}
\keyword{crop}
\keyword{rsf,}
