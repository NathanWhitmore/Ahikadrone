\name{unsupervised_pca}
\alias{unsupervised_pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unsupervised vegetation classification using principle component analysis
}
\description{
This function undertakes unsupervised vegetation classification using a principle component analysis.  The function will automatically request the user to (1) identify an exterior folder containing the pre-processed imagery, and (2) the folder which will contain the function's outputs. The function will report messages indicating progress through 9 steps:

\itemize{
\item \code{Step 1}: raster stacking starting
\item \code{Step 2}: raster shrinking starting
\item \code{Step 3}: PCA
\item \code{Step 4}: choose dimensions based on screeplot
\item \code{Step 5}: hierarchical clustering to guide number of clusters
\item \code{Step 6}: choose number of clusters
\item \code{Step 7}: random forest extrapolation of clustering
\item \code{Step 8}: de-noising using modal values
\item \code{Step 9}: polygon writing
}

At Step 4 and 6 the user will be asked for input. At Step 4 they will need to examine the screeplot and determine the number of dimensions to include. This can be done by examining the plot and deciding which is the last informative dimension to include. At Step 6 the user will be asked for the number of clusters that the dendrogram will be divided into. The user should make this assessment by identifying the major branches present and at which point further division becomes unwanted.

}


\usage{
unsupervised_pca(shrink = 10, nsamples = 5000, denoise = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{



\item{shrink:}{
A number representing how pixels in each direction (horizontally and vertically) should be aggregated together. A shrink of 10 will result in a raster image 100 x smaller than the original. The default is set to 10. Failure to shrink may result in the computer crashing as it may lack the RAM necessary for the calculations associated with very large rasters.
}

\item{denoise:}{
A number representing a window of pixels in each direction (horizontally and vertically) that will be used to smooth the local neighbourhood using modal values. This will effectively lessen the noise present in the raster. A default of 5 is used.
}
\item{nsamples:}{
The number of samples used to build the PCA through subsampling. A default of 5000 is used. Higher values will increase processing time.
}

}

\details{
The function uses PCA based on a sample defined by \code{nsamples}. The number of meaningful PCA dimensions to include are defined by the user based on a scree plot produced by the function. Hierarchical clustering is then used to visualise the likely hierarchy of clusters (through a dendrogram). Again the user is asked to determine the number of clusters, these will be shown on the dendrogram. The PCA algorithm is then extended to the entire raster using a specialised PCA algorithm (\code{rasterPCA}), while the clustering algorithm is then estimated using a random forests algorithm (as a euclidean distance matrix will exceed the available RAM).


One automatically named shape file  will be returned: the name will also record the number of dimensions and clusters used in its manufacture.
}

\examples{
unsupervised_pca(shrink = 10, nsamples = 5000, denoise = 5)
}
