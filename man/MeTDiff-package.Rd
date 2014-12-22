\name{MeTDiff-package}
\alias{MeTDiff-package}
\alias{MeTDiff}
\docType{package}
\title{
MeTDiff
}
\description{
The package is developed for the differential analysis for MeRIP-seq data of two experimental conditions to unveil the dynamics in post-transcriptional regulation of the RNA methylome. The MeTDiff R-package explicitly models the reads variation in data and also devices a more power likelihood ratio test for differential methylation site prediction. It accepts and statistically supports multiple biological replicates, internally removes PCR artifacts and multi-mapping reads, outputs exome-based binding sites (RNA methylation sites) and detects differential post-transcriptional RNA modification sites between two experimental conditions in term of percentage rather the absolute amount. The package is still under active development, and we welcome all biology and computation scientist for all kinds of collaborations and communications. Please feel free to contact Dr. Xiaodong Cui <xiaodong.cui@outlook.com> if you have any questions.
}
\details{
\tabular{ll}{
Package: \tab MeTDiff\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-18\cr
License: \tab GPL-2 \cr
}

}

\references{
Xiaodong Cui, Jia Meng, Manjeet K. Rao, Yidong Chen and Yufei Huang. "MeTDiff: a Novel Differential RNA Methylation Analysis for MeRIP-Seq Data."
}

\keyword{ Statistical Inference }

\examples{
# For usage, please check the main function with: 
?metdiff
}
