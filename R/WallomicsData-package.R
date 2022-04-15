#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#' Ecotype
#' @aliases Ecotype
#' @docType data
#' @description
#' The \strong{Ecotype factor} identifies the genotype specifically designed for a
#' given ecosystem of the \emph{A. thaliana} from which the studied sample comes from. We have a population of reference as well as
#' 4 newly-described Pyrenean populations, namely:
#' \itemize{
#' \item Columbia, denoted \emph{Col} (originating from Poland, acts as the reference ecotype)
#' \item Grip, denoted \emph{Grip}
#' \item Herran, denoted \emph{Hern}
#' \item L’Hospitalet-près-l’Andorre, denoted \emph{Hosp}
#' \item Chapelle Saint Roch, denoted \emph{Roch}
#' }
#' @usage data("Ecotype")
#' @format A factor with 5 levels of \emph{A. thaliana} genotypes.
#' @examples
#' # Load the data
#' data("Ecotype")
#'
#' # Count how many samples are in each group
#' table(Ecotype)
#'
#' @source \doi{10.3390/cells9102249}
"Ecotype"

#' Temperature
#' @aliases Temperature
#' @docType data
#' @description
#' The \strong{Temperature factor} identifies the temperature at which the studied sample was exposed
#' all along its growth, either 22°C (optimal condition) or 15°C (high altitude condition).
#' @usage data("Temperature")
#' @format A factor with 2 levels.
#' @examples
#' # Load the data
#' data("Temperature")
#'
#' # Count how many samples are in each group
#' table(Temperature)
#'
#' @source \doi{10.3390/cells9102249}
"Temperature"

#' Altitude Cluster
#' @docType data
#' @description
#' The \strong{Altitude Cluster factor} identifies the environment height
#' from which is originated a given plant from the sample under study, either high altitude (denoted \emph{High}),
#' moderate altitude (\emph{Low}) or the reference group's environment height (\emph{Col}, the lowest of all).
#' @usage data("Altitude_Cluster")
#' @format A factor with 3 levels.
#' @examples
#' # Load the data
#' data("Altitude_Cluster")
#'
#' # Count how many samples are in each group
#' table(Altitude_Cluster)
#'
#' @source \doi{10.3390/cells9102249}
"Altitude_Cluster"

#' Genetic Cluster
#' @docType data
#' @description
#' The \strong{Genetic Cluster factor} identifies the genetic group
#' from which comes from the studied sample, either \strong{Genetics Cluster 1} (constitued of \emph{Grip} and \emph{Roch} genotypes),
#' \strong{Genetics Cluster 2} (\emph{Hern} and \emph{Hosp} genotypes) or \strong{Genetics Cluster 3} (\emph{Col} genotype).
#' See \link{Ecotype} for more information on genotypes.
#' @usage data("Genetic_Cluster")
#' @format A factor with 3 levels.
#' @examples
#' # Load the data
#' data("Genetic_Cluster")
#'
#' # Count how many samples are in each group
#' table(Genetic_Cluster)
#'
#' @source \doi{10.3390/cells9102249}
"Genetic_Cluster"

#' Phenomics Rosettes
#' @docType data
#' @description
#' A dataset containing phenotypic variables measured on rosettes of
#' five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Phenomics_Rosettes")
#' @format A data frame with 30 rows and 5 variables:
#' \itemize{
#'   \item \strong{Mass}: rosette mass (g)
#'   \item \strong{Diameter}: rosette diameter (cm)
#'   \item \strong{Leaves_number}: total number of leaves
#'   \item \strong{Density}: rosette density (g/cm²)
#'   \item \strong{Area}: projected rosette area (cm²)
#' }
#' @examples
#' # Load the data
#' data("Phenomics_Rosettes")
#'
#' # Look at simple statistics
#' dim(Phenomics_Rosettes)
#' summary(Phenomics_Rosettes)
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # A graphical representation: Leaves number distribution
#' plot(x = as.factor(substr(row.names(Phenomics_Rosettes), 1, 7)),
#'      y = Phenomics_Rosettes$Leaves_number, col = "white", lty = 0,
#'      xlab = "Genotype x Temperature groups",
#'      ylab = "Number of rosette leaves",
#'      main = "Rosette leaves' distribution by genotype and growth temperature"
#'      )
#' grid()
#' points(x = as.factor(substr(row.names(Phenomics_Rosettes), 1, 7)),
#'        y = Phenomics_Rosettes$Leaves_number, type = "p", pch = 19, lwd = 5,
#'        col = colors)
#'
#' @source \doi{10.3390/cells9102249}
"Phenomics_Rosettes"

#' Phenomics Stems
#' @docType data
#' @description
#' A dataset containing phenotypic variables measured on floral stems of
#' five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Phenomics_Stems")
#' @format A data frame with 30 rows and 4 variables:
#' \itemize{
#'   \item \strong{Mass}: floral stems mass (g)
#'   \item \strong{Diameter}: floral stems diameter (mm)
#'   \item \strong{Length}: length of the floral stems (cm)
#'   \item \strong{Number_lateral_stems}: number of lateral stems)
#' }
#' @examples
#' # Load the data
#' data("Phenomics_Stems")
#'
#' # Look at simple statistics
#' dim(Phenomics_Stems)
#' summary(Phenomics_Stems)
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # A graphical representation: Lateral stems distribution
#' plot(x = as.factor(substr(row.names(Phenomics_Stems), 1, 7)),
#'      y = Phenomics_Stems$Number_lateral_stems, col = "white", lty = 0,
#'      xlab = "Genotype x Temperature groups",
#'      ylab = "Number of lateral stems",
#'      main = "Lateral stems' distribution by genotype and growth temperature"
#'      )
#' grid()
#' points(x = as.factor(substr(row.names(Phenomics_Stems), 1, 7)),
#'        y = Phenomics_Stems$Number_lateral_stems, type = "p", pch = 19, lwd = 5,
#'        col = colors)
#'
#' @source \doi{10.3390/cells9102249}
"Phenomics_Stems"

#' Metabolomics Rosettes
#' @docType data
#' @description
#' A dataset containing metabolomics variables measured on rosettes of
#' five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Metabolomics_Rosettes")
#' @format A data frame with 30 rows and 6 variables:
#' \itemize{
#'   \item \strong{Pectin_RGI}: Rhamnogalacturonan I (µg/100mg)
#'   \item \strong{Pectin_HG}: Homogalacturonan (µg/100mg)
#'   \item \strong{XG}: Xyloglucan (µg/100mg)
#'   \item \strong{Pectin_linearity}: Linearity of pectin (Ratio)
#'   \item \strong{Contribution_RG}: Contribution of rhamnogalacturonan to pectin population (Ratio)
#'   \item \strong{RGI_branching}: Branching of Rhamnogalacturonan I (Ratio)
#' }
#' @examples
#' # Load the dataset
#' data("Metabolomics_Rosettes")
#'
#' # Look at simple statistics
#' summary(Metabolomics_Rosettes)
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # A graphical representation
#' plot(x = as.factor(substr(row.names(Metabolomics_Rosettes), 1, 7)),
#'      y = Metabolomics_Rosettes$Pectin_linearity, col = "white", lty = 0,
#'      xlab = "Genotype x Temperature groups",
#'      ylab = "Pectin linearity (Ratio)",
#'      main = "Pectin linearity distribution by genotype and growth temperature")
#' grid()
#' abline(h = 1, lty = 2)
#' points(x = as.factor(substr(row.names(Metabolomics_Rosettes), 1, 7)),
#'        y = Metabolomics_Rosettes$Pectin_linearity, type = "p", pch = 19, lwd = 5,
#'        col = colors)
#'
#' @source \doi{10.3390/cells9102249}
"Metabolomics_Rosettes"

#' Metabolomics Stems
#' @docType data
#' @description
#' A dataset containing metabolomics variables measured on floral stems of
#' five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Metabolomics_Stems")
#' @format A data frame with 30 rows and 6 variables:
#' \itemize{
#'   \item \strong{Pectin_RGI}: Rhamnogalacturonan I (µg/100mg)
#'   \item \strong{Pectin_HG}: Homogalacturonan (µg/100mg)
#'   \item \strong{XG}: Xyloglucan (µg/100mg)
#'   \item \strong{Pectin_linearity}: Linearity of pectin (Ratio)
#'   \item \strong{Contribution_RG}: Contribution of rhamnogalacturonan to pectin population (Ratio)
#'   \item \strong{RGI_branching}: Branching of Rhamnogalacturonan I (Ratio)
#' }
#' @examples
#' # Load the dataset
#' data("Metabolomics_Stems")
#'
#' # Look at simple statistics
#' summary(Metabolomics_Stems)
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # A graphical representation
#' plot(x = as.factor(substr(row.names(Metabolomics_Stems), 1, 7)),
#'      y = Metabolomics_Stems$Pectin_linearity, col = "white", lty = 0,
#'      xlab = "Genotype x Temperature groups",
#'      ylab = "Pectin linearity (Ratio)",
#'      main = "Pectin linearity distribution by genotype and growth temperature")
#' grid()
#' abline(h = 1, lty = 2)
#' points(x = as.factor(substr(row.names(Metabolomics_Stems), 1, 7)),
#'        y = Metabolomics_Stems$Pectin_linearity, type = "p", pch = 19, lwd = 5,
#'        col = colors)
#'
#' @source \doi{10.3390/cells9102249}
"Metabolomics_Stems"

#' Proteomics Rosettes Cell Wall
#' @docType data
#' @description
#' A dataset containing the identification and quantification of Cell Wall Proteins (CWPs) performed using LC-MS/MS
#' analysis on rosettes of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for additional information.
#' @usage data("Proteomics_Rosettes_CW")
#' @format A data frame with 30 rows and 364 variables.
#' @examples
#' # Load the dataset
#' data("Proteomics_Rosettes_CW")
#'
#' # Look at data frame dimensions
#' dim(Proteomics_Rosettes_CW)
#'
#' # Look at the first rows and columns
#' head(Proteomics_Rosettes_CW[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on proteomics
#' res.pca <- prcomp(Proteomics_Rosettes_CW, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Rosettes Cell Wall Proteomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Proteomics_Rosettes_CW"

#' Proteomics Stems Cell Wall
#' @docType data
#' @description
#' A dataset containing the identification and quantification of Cell Wall Proteins (CWPs) performed using LC-MS/MS
#' analysis on floral stems of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for additionnal information.
#' @usage data("Proteomics_Stems_CW")
#' @format A data frame with 30 rows and 414 variables.
#' @examples
#' # Load the dataset
#' data("Proteomics_Stems_CW")
#'
#' # Look at data frame dimensions
#' dim(Proteomics_Stems_CW)
#'
#' # Look at the first rows and columns
#' head(Proteomics_Stems_CW[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on proteomics
#' res.pca <- prcomp(Proteomics_Stems_CW, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Stems Cell Wall Proteomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Proteomics_Stems_CW"

#' Transcriptomics Rosettes Cell Wall
#' @docType data
#' @description
#' A dataset containing the transcripts encoding Cell Wall Proteins (CWPs)
#' sorted from the 19 763 transcripts (see \link{Transcriptomics_Rosettes}) obtained by RNA-seq performed,
#' according to the standard Illumina protocols, on rosettes of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Transcriptomics_Rosettes_CW")
#' @format A data frame with 30 rows and 364 variables.
#' @examples
#' # Load the dataset
#' data("Transcriptomics_Rosettes_CW")
#'
#' # Look at data frame dimensions
#' dim(Transcriptomics_Rosettes_CW)
#'
#' # Look at the first rows and columns
#' head(Transcriptomics_Rosettes_CW[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on transcriptomics
#' res.pca <- prcomp(Transcriptomics_Rosettes_CW, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Rosettes Cell Wall Transcriptomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Transcriptomics_Rosettes_CW"

#' Transcriptomics Stems Cell Wall
#' @docType data
#' @description
#' A dataset containing the transcripts encoding Cell Wall Proteins (CWPs) sorted
#' from the 22 570 transcripts (see \link{Transcriptomics_Stems}) obtained by RNA-seq performed,
#' according to the standard Illumina protocols,
#' on floral stems of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Transcriptomics_Stems_CW")
#' @format A data frame with 30 rows and 414 variables.
#' @examples
#' # Load the dataset
#' data("Transcriptomics_Stems_CW")
#'
#' # Look at data frame dimensions
#' dim(Transcriptomics_Stems_CW)
#'
#' # Look at the first rows and columns
#' head(Transcriptomics_Stems_CW[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on transcriptomics
#' res.pca <- prcomp(Transcriptomics_Stems_CW, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Stems Cell Wall Transcriptomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Transcriptomics_Stems_CW"

#' Transcriptomics Rosettes
#' @aliases Transcriptomics_Rosettes
#' @docType data
#' @description
#' A dataset containing all the transcripts obtained by RNA-seq performed, according to the standard Illumina protocols,
#' on rosettes of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Transcriptomics_Rosettes")
#' @format A data frame with 30 rows and 19763 variables.
#' @examples
#' # Load the dataset
#' data("Transcriptomics_Rosettes")
#'
#' # Look at data frame dimensions
#' dim(Transcriptomics_Rosettes)
#'
#' # Look at the first rows and columns
#' head(Transcriptomics_Rosettes[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on transcriptomics
#' res.pca <- prcomp(Transcriptomics_Rosettes, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Rosettes' Transcriptomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Transcriptomics_Rosettes"

#' Transcriptomics Stems
#' @aliases Transcriptomics_Stems
#' @docType data
#' @description
#' A dataset containing all the transcripts obtained by RNA-seq performed, according to the standard Illumina protocols,
#' on floral stems of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for more information.
#' @usage data("Transcriptomics_Stems")
#' @format A data frame with 30 rows and 22570 variables.
#' @examples
#' # Load the dataset
#' data("Transcriptomics_Stems")
#'
#' # Look at data frame dimensions
#' dim(Transcriptomics_Stems)
#'
#' # Look at the first rows and columns
#' head(Transcriptomics_Stems[,c(1:10)])
#'
#' # Create a colors' vector
#' colors <- c(rep("#A6CEE3",3), rep("#1F78B4",3), rep("#B2DF8A",3), rep("#33A02C",3),
#'             rep("#FB9A99",3), rep("#E31A1C",3), rep("#FDBF6F",3), rep("#FF7F00",3),
#'             rep("#CAB2D6",3), rep("#6A3D9A",3))
#'
#' # PCA on transcriptomics
#' res.pca <- prcomp(Transcriptomics_Stems, center = TRUE, scale. = TRUE)
#' plot(res.pca$x[,"PC1"], res.pca$x[,"PC2"], pch = 19, xlab = "PC1", ylab = "PC2", lwd = 5,
#'      main = "Individuals' plot (1 x 2) - PCA on Stems' Transcriptomics",
#'      col = colors)
#' text(res.pca$x[,"PC1"], res.pca$x[,"PC2"], labels = row.names(res.pca$x), cex = 0.8, pos = 3)
#'
#' @source \doi{10.3390/cells9102249}
"Transcriptomics_Stems"

#' Metadata
#' @docType data
#' @description
#' Bioinformatics Annotation and description, using the WallProtDB database, of all the Cell Wall Proteins (CWPs)
#' identified on rosettes and floral stems of five \emph{A. thaliana} genotypes at two growth temperatures.
#' See \link{Ecotype} and \link{Temperature} for additionnal information.
#' @usage data("Metadata")
#' @format A data frame with 474 rows and 4 variables:
#' \itemize{
#'   \item \strong{Acc_number}: GenBank accession number (gene name)
#'   \item \strong{Functional_classes}: Functional classes of the CWPs
#'   \item \strong{Protein_families}: Protein families of the CWPs
#'   \item \strong{Putative_functions}: Putative functions of the CWPs
#' }
#' @examples
#' # Load the dataset
#' data("Metadata")
#'
#' # Look at the dataset's dimensions
#' dim(Metadata)
#' head(Metadata)
#'
#' # How many functional classes ?
#' table(Metadata$Functional_classes)
#'
#' # How many protein families ?
#' table(Metadata$Protein_families)
#'
#' @source \doi{10.3390/cells9102249}
"Metadata"
