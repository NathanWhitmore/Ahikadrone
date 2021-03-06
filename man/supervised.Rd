\name{supervised}
\alias{supervised}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Undertaking supervised (predictive) vegetation classification using random forests
}
\description{
This function undertakes supervised vegetation classification from DJI multispectral drone imagery that has been already pre-processed through the DJI propriety software. In addition to the two arguments (shrink and seed), the function will automatically request the user to (1) identify exterior folder containing the pre-processed imagery, (2) file the ESRI shapefile that will be used for training, and (3) the folder which will contain the function`s outputs.

The function will report messages indicating progress through 6 steps:

\itemize{
\item \code{Step 1}: raster stacking starting
\item \code{Step 2}: raster shrinking starting
\item \code{Step 3}: extracting training layers
\item \code{Step 4}: undertaking PCA
\item \code{Step 5}: undertaking random forests
\item \code{Step 6}: writing polygons
}

\bold{Note:} Ideally each training class should contain a minimum of 200 pixels across all polygons for each vegetation type. For example, if a 20x shrink is planned, a minimum of 4000 pixels would be required to be preidentified in order to recognise a specific vegetation class e.g. "pine". This is because, in order to prevent class imbalances from occurring, the algorithm needs to sample the training dataset by randomly sampling 200 pixels within each training class after shrinking has occurred. Failure to do this will result in a class imbalance. For rare vegetation classes this may be unavoidable.
}

\usage{
supervised(shrink = 10)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{
a number representing how pixels in each direction (horizontally and vertically) should be aggregated together. A shrink of 10 will result in a raster image 100 x smaller than the original. The default is set to 10. Failure to shrink may result in the computer crashing as it may lack the RAM necessary for the calculations associated with very large rasters.
}

}
\details{
The function operates by first aggregating the raster images (using median values) to a smaller size so that they are manageable by a desktop PC, then reduces all multispectral channels (blue, green, red, red edge, and near infra-red) and derived index layers (NDVI, GNDVI, LCI, NDRE, and OSAVI) into uncorrelated explanatory variables using principle component analysis (note: visible light wavelengths remain important for discrimination of non-photosynthetic land-use classes). These explanatory variables are then used to predict the response variable (vegetation type) given by the training shape file using the random forest machine learning algorithm and tuned using 5-fold cross validation.
Two outputs are produced (1) a text summary of the accuracy of the random forest algorithm and identification of any vegetation classes under-represented in the training shape file and (2) the ESRI shape files associated with the vegetation classes. Any polygons 1 square metre or less are purposely removed to prevent the shapefiles by being overwhelmed by noise. Based on the class imbalances identified in the analysis the user may want to increase the size of the polygons in the training shape files if they appear to be under-represented. Ideally, the training file should identify examples of vegetation classes that are both representative of the variation present in the image and relatively homogeneous in terms of their composition (e.g. polygons identifying pine should avoid including another species such as manuka).
}


\examples{
supervised(shrink = 10)
}

