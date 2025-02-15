% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WallomicsData-package.R
\docType{data}
\name{Phenomics_Rosettes}
\alias{Phenomics_Rosettes}
\title{Phenomics Rosettes}
\format{
A data frame with 30 rows and 5 variables:
\itemize{
\item \strong{Mass}: rosette mass (g)
\item \strong{Diameter}: rosette diameter (cm)
\item \strong{Leaves_number}: total number of leaves
\item \strong{Density}: rosette density (g/cm²)
\item \strong{Area}: projected rosette area (cm²)
}
}
\source{
\doi{10.3390/cells9102249}
}
\usage{
data("Phenomics_Rosettes")
}
\description{
A dataset containing phenotypic variables measured on rosettes of
five \emph{A. thaliana} genotypes at two growth temperatures.
See \link{Ecotype} and \link{Temperature} for more information.
}
\examples{
# Load the data
data("Phenomics_Rosettes")

# Look at simple statistics
dim(Phenomics_Rosettes)
summary(Phenomics_Rosettes)

# Create a colors' vector
colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
            rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
            rep("#CAB2D6",3), rep("#6A3D9A",3))

# A graphical representation: Leaves number distribution
plot(x = as.factor(substr(row.names(Phenomics_Rosettes), 1, 7)),
     y = Phenomics_Rosettes$Leaves_number, col = "white", lty = 0,
     xlab = "Genotype x Temperature groups",
     ylab = "Number of rosette leaves",
     main = "Rosette leaves' distribution by genotype and growth temperature"
     )
grid()
points(x = as.factor(substr(row.names(Phenomics_Rosettes), 1, 7)),
       y = Phenomics_Rosettes$Leaves_number, type = "p", pch = 19, lwd = 5,
       col = colors)

}
\keyword{datasets}
