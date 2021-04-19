\name{unsupervised}
\alias{unsupervised}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unsupervised vegetation classification using affinity propagation and k-means clustering
}
\description{
This function undertakes unsupervised vegetation classification using affinity propagation and k-means clustering.  The function will automatically request the user to (1) identify exterior folder containing the pre-processed imagery, and (2) the folder which will contain the function's outputs. The function will report messages indicating progress through 6 steps:

\itemize{
\item \code{Step 1}: raster stacking starting
\item \code{Step 2}: raster shrinking starting
\item \code{Step 3}: image segmentation starting
\item \code{Step 4}: affinity propagation starting (clustering based on physical proximity)
\item \code{Step 5}: k-means clustering starting (clustering based on shared attributes)
\item \code{Step 6}: polygon writing
}

These messages are important because it will enable the user understand if a failure of the function occurs because the raster has not be shrunk sufficiently.

}
\usage{
unsupervised(type, super, clusters, shrink)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{type}{
One of 7 analysis types:

\itemize{
\item \code{visible}: the 3-channel visible spectrum (blue, green, red).
\item \code{false_colour_nir}: a 3-channel red-shifted spectrum (green, red, near infra-red) useful for vegetation discrimination.
\item \code{ndvi}: Normalized Difference Vegetation Index. The most common index used in agriculture. It characterises the density of vegetation and allows farmers to assess germination, growth, the presence of weeds or diseases, as well as to predict the productivity of the fields.
\item \code{gdvi}: Green Normalized Difference Vegetation Index. It is similar to NDVI except that instead of the red spectrum it measures the green spectrum in the range from 0.54 to 0.57 microns. This is an indicator of the photosynthetic activity of the vegetation cover; it is most often used in assessing the moisture content and nitrogen concentration in plant leaves according to multispectral data which do not have an extreme red channel. Compared to the NDVI index, it is more sensitive to chlorophyll concentration. It is commonly used in assessing depressed and aged vegetation.
\item \code{lci}: Leaf Chlorophyll Index. It is used to assess chlorophyll content in areas of complete leaf coverage.
\item \code{ndre}: Normalized Difference Red Edge Index. This index is sensitive to chlorophyll content in leaves against soil background effects. This index can only be formulated when the red edge band is available
\item \code{osavi}: Optimized Soil-Adjusted Vegetation Index. This index maps variability in canopy density. It is not sensitive to soil brightness (when different soil types are present.
}

}
  \item{super}{
The approximate number of segments (superpixels) the visible image will be divided into. A default of 200 is recommended. However a range of 80 -- 300 may be sensible. High values will increase processing time.
}
  \item{clusters}{
The number of land use clusters (vegetation types) you are wanting to discriminate. If the number of specified clusters exceeds the number of clusters returned by the affinity propagation algorithm, the maximum number of clusters allowable by the affinity propagation algorithm will be used.
}
  \item{shrink}{
A number representing how pixels in each direction (horizontally and vertically) should be aggregated together. A shrink of 10 will result in a raster image 100 x smaller than the original. The default is set to 10. Failure to shrink may result in the computer crashing as it may lack the RAM necessary for the calculations associated with very large rasters.
}
}
\details{
Two automatically named polygon files will be returned: one containing the clustered polygons, and the other containing the polygons of the superpixels (image segments).
An automatically named elbow plot will be produced as a .png file. While there is no set rule for determining the number of clusters, many people look for an inflection point (the elbow, or knee) in this graph as a guide to determining the correct number of clusters.
}


\examples{
unsupervised(type = "ndre",
            clusters = 6,
            super = 200,
            shrink = 10)
  }


