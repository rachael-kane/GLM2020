% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cofactor.pca.cor.R
\name{cofactor.pca.cor}
\alias{cofactor.pca.cor}
\title{Correlation between cofactors and principal components.}
\usage{
cofactor.pca.cor(U, G)
}
\arguments{
\item{U}{A numeric matrix containing user-specified cofactors. Dimensions are n rows (individuals) by t columns (cofactors).}

\item{G}{A numeric matrix containing genotype data. Dimensions are n rows (individuals) by m columns (genetic markers).}
}
\value{
A list of 1 or 3 objects.

U unspecified: 1 object.
 $cov, a numeric matrix containing all principal components and individual scores.

U specified: 3 objects.
 $orig_pc, a numeric matrix containing all original principal components
 $cov, a numeric matrix containing user-specified cofactors and retained principal components.
 $removed, a matrix indicating which principal components were removed.
}
\description{
Test for correlations between user-specified cofactors and principal components calculated from genotype data. Automatically remove principal components linearly dependent (correlated) with user-specified cofactors.
}
\details{
When U is unspecified, cofactor.pca.cor will return a list of 1 object.
With U unspecified, function will carry out principal components analysis identically to the native R function prcomp(),
 and cofactor.pca.cor will return principal components scores in $cov.
$cov is a numeric matrix containing all principal components and individual scores.
 Dimensions are n rows (individuals) by t columns (principal components).

When U is specified, cofactor.pca.cor will return a list of 3 objects.
$orig_pc is a numeric matrix containing all original principal components and individual scores.
$cov is a numeric matrix containing user-specified cofactors and all principal components not correlated with the
 user-specified cofactors. Dimensions are n rows (individuals) by t columns (cofactors).
$removed is a character matrix indicating which principal components were removed.

The $cov matrix is intended for use as the "C" argument in the GWASbyGLM function included in this package.

Type vignette("GLM2020_tutorial") for example use.
}
