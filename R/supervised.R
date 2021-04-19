


supervised <- function(seed = TRUE, shrink = 10){

  if(seed == TRUE) {set.seed(27)} else {set.seed(NULL)}


  # find folder
  dir <- choose.dir(default = "", caption = "Select folder which contains DJI drone imagery\nNote: only select the exterior folder")

  # set path
  path <- paste0(dir, "\\map\\")

  # output name
  output.name <- sub(".*\\\\", "", dir)

  # finding training file
  training.file <- choose.files(default = "", caption = "Select the ESRI shapefile (.shp) that will be used for training")

  # output folder
  out.dir <- choose.dir(default = "", caption = "Select folder where the output files are to be put")

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

  # make 3-d brick all (3 minutes)
  region.brick <- brick(blue, green, red, red_edge, nir, ndvi, gndvi, lci, ndre, osavi)

  message("Step 2: raster shrinking starting")

  # make lower resolution raster (10 x smaller) 12 minutes with fun = median
  region.brick <- aggregate(region.brick, fact = shrink, fun = median)

  # read files
  supervised <- st_read(training.file)

  # add transform into longlat
  supervised <- st_transform(supervised, crs = "+proj=longlat")

  # add id
  supervised$id <- 1:nrow(supervised)

  # make count (IMPORTANT)
  supervised$code <- as.numeric(as.factor(supervised$VegType))

  # unique vege names and number
  vege.names.no <- length(unique(supervised$code))
  vege.types <- sort(unique(supervised$VegType))

  # change to spatial type
  supervised <- as_Spatial(supervised)

  # mask
  masked.brick <- mask(region.brick, supervised, sp = TRUE)

  print(vege.names.no)
  print(vege.types)
  class(vege.types)


  # for loop to extract as separate layers
  ground_truth <- NULL
  my.extract <- NULL

  message("Step 3: extracting training layers")

  for (i in 1:vege.names.no)  {
    # get correct layer
    my.mask <- supervised[supervised$code == i,]
    masked <- raster::mask(region.brick, my.mask)

    # get sample points (fixed at 100)
    ground_truth <- as.data.frame(sampleRandom(masked, 200, xy= TRUE))
    ground_truth$vege <- as.character(vege.types[i])

    # rbind to original
    my.extract <- rbind(ground_truth , my.extract)
  }

  # change to factor vege.types
  levels(my.extract$vege) <- factor(vege.types)

  message("Step 4: undertaking PCA")

  # run pca to get non-correlated variables
  pca <- princomp(my.extract[,-c(1,2,length(my.extract))])

  # pca data frame
  pca.df <- cbind(pca$scores, data.frame(vege_type = my.extract[,length(my.extract)]))

  message("Step 5: undertaking random forests")

  # randomforest
  cross.val <- trainControl(method = "cv", number = 5)
  train.rf <- train(vege_type ~ .,
                    data = pca.df,
                    method = "rf",
                    trControl = cross.val)

  # confusion matrix
  cm <- confusionMatrix(train.rf)

  # confusion matrix output
  df.con <- as.data.frame(cm$table)
  wide.con <- df.con %>% pivot_wider(names_from = Reference, values_from = Freq)
  diag <- diag(as.matrix(wide.con[,-1]))

  wide.con <- wide.con %>% adorn_totals()

  last.line <- tail(wide.con,1)
  last.line <- last.line[,- 1]

  perfect <- diag/ last.line

  # make list
  output <- list()

  output$random.forest <- train.rf
  output$confusion.matrix <- cm
  output$perfect <- perfect
  output$imbalance <- last.line /max(last.line)

  print(output)

  # save output
  capture.output(output, file = paste0(out.dir, "\\", output.name,"_supervised", ".txt"), append = FALSE)


  # make predictions based on pca
  pca.df$test <- predict(train.rf,  pca.df)


  # PCA on entire raster
  pca.raster <- rasterPCA(region.brick, nSamples = 1000)

  pca.raster.df <- raster::as.data.frame(pca.raster$map, xy = TRUE)
  pca.raster.df <- na.omit(pca.raster.df)

  # rename
  colnames(pca.raster.df) <- c("x", "y", paste0("Comp.", 1:ncol(pca$scores)))

  message("prediction")

  # prediction
  pred.values <- raster::predict(train.rf,  pca.raster.df)

  # make new df
  my.data <- pca.raster.df

  print(head(my.data))

  # layer names
  my.data$vege_pred <- pred.values

  # change to raster

  # change data to numeric
  my.data$vege <- as.numeric(my.data$vege_pred)

  print(head(my.data))

  message("change df to raster")

  # change df to raster
  my.data <- my.data[, c(1,2,length(my.data))]

  # make into raster
  dfr <- rasterFromXYZ(my.data)  #Convert first two columns as lon-lat and third as value
  crs(dfr) <- crs(region.brick)

  print(dfr)

  message("change dfr")
  # make into factor
  dfr <- raster::as.factor(dfr)

  print(dfr)

  message("make table to describe factors A")


  x <- raster::as.data.frame(levels(dfr))

  print(head(x))

  message("make table to describe factors B")

  x$vege_type <- vege.types

  message("make table to describe factors C")

  raster::levels(dfr) <- x

  message("Step 6: writing polygons")

  # change to polygon
  orig.poly <- rasterToPolygons(dfr, na.rm=TRUE, dissolve=TRUE)

  decrumbed <- drop_crumbs(orig.poly, set_units(1, m^2))

  # change to sf for easier handling
  vege.sf <- st_as_sf(decrumbed)

  # add actual levels to sf as these have been lost
  vege.sf$vege_type <- vege.types

  # write to shp for QGIS
  st_write(vege.sf, paste0(out.dir, "\\", output.name,"_supervised",".shp"), driver = "ESRI Shapefile", append=FALSE)

}
