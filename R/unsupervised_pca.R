
unsupervised_pca <- function(shrink = 10, nsamples = 5000, denoise = 5){

  if(denoise %% 2 == 0) stop('Denoise is not odd')


  # find folder
  dir <- choose.dir(default = "", caption = "Select folder which contains DJI drone imagery\nNote: only select the exterior folder")

  # output name
  output.name <- sub(".*\\\\", "", dir)

  # output folder
  out.dir <- choose.dir(default = "", caption = "Select folder where the output files are to be put")

  # set path
  path <- paste0(dir, "\\map\\")

  # spectral
  blue <- raster(paste0(path , "result_Blue.tif"))
  green <- raster(paste0(path , "result_Green.tif"))
  red <- raster(paste0(path , "result_Red.tif"))
  red_edge <- raster(paste0(path , "result_RedEdge.tif"))
  nir <- raster(paste0(path , "result_NIR.tif"))

  # derived
  ndvi <- raster(paste0(path , "\\index_map\\", "NDVI.tif"))
  gndvi <- raster(paste0(path , "\\index_map\\", "GNDVI.tif"))
  lci <- raster(paste0(path , "\\index_map\\", "LCI.tif"))
  ndre <- raster(paste0(path , "\\index_map\\", "NDRE.tif"))
  osavi <- raster(paste0(path , "\\index_map\\", "OSAVI.tif"))

  message("Step 1: raster stacking starting")

  region.brick <- brick(blue, green, red, red_edge, nir, ndvi, gndvi, lci, ndre, osavi)

  message("Step 2: raster shrinking starting")

  # make lower resolution raster (10 x smaller) 12 minutes with fun = median
  region.brick <- aggregate(region.brick, fact = shrink, fun = median)

  message("Step 3: PCA")


  # get sample points
  df_sample <- as.data.frame(sampleRandom(region.brick, nsamples, xy= TRUE))
  pca_sample <-  princomp(df_sample[,-c(1,2)])
  pca.sample.df <- as.data.frame(cbind(pca_sample$scores))

  message("Step 4: choose dimensions based on screeplot")

  print(fviz_eig(pca_sample))

  my.dim <- as.numeric(readline(prompt="Enter number of dimensions: "))

  # dimension reduction based on user
  pca.sample.df <- pca.sample.df[, 1: my.dim]

  message("Step 5: hierarchical clustering to guide number of clusters")

  # hierarchical clustering to guess number of clusters
  # cult.results <- HCPC(pca.sample.df, min = min.clusters , max = max.clusters , graph = FALSE)


#  cult.results <- HCPC(pca.sample.df, min = min.clusters , max = max.clusters , graph = FALSE)
  #  data.clusters$clust <- factor(data.clusters$clust)
  #  k.clust <- length(unique(data.clusters$clust))

  dendro <- pca.sample.df %>% dist(method = "euclidean") %>%
    hclust(method = "ward.D")

  # initial plot
  plot(dendro, labels = FALSE)

  message("Step 6: choose number of clusters")
  my.cuts <- as.numeric(readline(prompt="Enter number of clusters: "))

  # replot
  plot(dendro, labels = FALSE)
  rect.hclust(dendro, k =  my.cuts)

  # cuts
  clust <- as.factor(cutree(dendro, k= my.cuts))

  # new data frame
  data.clusters <-  cbind(pca.sample.df, clust)

  message("Step 7: Random forest extrapolation of clustering")

  # randomforest
  cross.val <- trainControl(method = "cv", number = 5)
  train.rf <- train(clust ~ .,
                    data = data.clusters,
                    method = "rf",
                    trControl = cross.val)

  # confusion matrix
  cm <- confusionMatrix(train.rf)
  cm

  # PCA on entire raster
  pca.raster <- rasterPCA(region.brick, nSamples = nsamples)
  pca.raster.df <- raster::as.data.frame(pca.raster$map, xy = TRUE)


  # rename
  colnames(pca.raster.df) <- c("x", "y", paste0("Comp.", 1:ncol(pca_sample$scores)))


  # prediction including handling of NAs
  na.free <- na.omit(pca.raster.df)
  pca.raster.df$clusters <- NA
  pca.raster.df$clusters[which(!is.na(pca.raster.df[,3]))] <-  predict(train.rf,  na.free)


  # make into raster
  my.clust <- setValues(region.brick[[1]], pca.raster.df$cluster)

  message("Step 8: De-noising using modal values")
  focused <- focal(my.clust, w=matrix(1/(denoise^2),nrow=denoise,ncol=denoise), fun = modal, pad = TRUE) * (denoise^2)

  message("Step 9: polygon writing")

  # change to polygon
  orig.poly <- st_as_stars(focused) %>% st_as_sf(merge=TRUE)

  # change to sf for easier handling
  pca.sf <- st_as_sf(orig.poly)
  pca.sf <- pca.sf %>% rename(pca.clusters = layer)


  # write to shp for QGIS
  suppressWarnings(st_write(pca.sf, paste0(out.dir, "\\", output.name, "_unsupervised_pca_", my.dim, "_dim_",
                                           my.cuts, "_clust",
                                           ".shp"),
                            driver = "ESRI Shapefile", append =FALSE))

  message("Process completed")

}
