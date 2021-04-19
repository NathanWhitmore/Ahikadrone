

unsupervised <- function(type, super = 200, clusters, shrink = 10){

  set.seed(27)


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

  if(type == "visible"){
    region.brick <- brick(blue, green, red)

  } else if (type == "false_colour_nir"){
    region.brick <- brick(green, red, nir)

  } else if (type == "ndvi"){
    region.brick <- brick(ndvi, ndvi, ndvi)

  } else if (type == "gndvi"){
    region.brick <- brick(gndvi, gndvi, gndvi)

  } else if (type == "lci"){
    region.brick <- brick(lci, lci, lci)

  } else if (type == "ndre"){
    region.brick <- brick(ndre, ndre, ndre)

  } else if (type == "osavi"){
    region.brick <- brick(osavi, osavi, osavi)

  } else {
    message(paste0("Type '", type, "' is not available as an option"))
  }

  message("Step 2: raster shrinking starting")

  # shrink resolution
  region.brick <- aggregate(region.brick, fact = shrink)

  # necessary for writing image to disk
  nrows <- region.brick@nrows
  ncols <- region.brick@ncols

  # write to disk
  long.name <- paste0(out.dir,"\\", output.name,"_", "unsupervised_", type, ".jpg")
  jpeg(long.name, width = ncols, height = nrows)
  plotRGB(region.brick, r = 3, g = 2, b = 1, stretch = "lin")
  dev.off()

  # read false colour
  my.image <- readImage(long.name)


  # super pixel calculator
  tot.cells <- nrows * ncols
  non.na <- region.brick[[1]]
  non.na[!is.na(non.na)] <- 1
  not.na <- freq(non.na)[1,2]
  conv.factor <- not.na/tot.cells


  message("Step 3: image segmentation starting")

  Region.slic <- superpixels(input_image = my.image,
                             method = "slic",
                             superpixel = super/conv.factor,
                             return_slic_data = TRUE,
                             return_labels = TRUE, write_slic = "",
                             verbose = FALSE)

  imageShow(Region.slic$slic_data)

  message("Step 4: affinity propagation starting")

  # actual segmentation using Affinity Propagation
  init <- Image_Segmentation$new()

  Region.seg <- init$spixel_segmentation(input_image = my.image,
                                         superpixel = super/conv.factor,
                                         AP_data = TRUE,
                                         use_median = TRUE,
                                         sim_wL = 3,
                                         sim_wA = 10,
                                         sim_wB = 10,
                                         sim_color_radius = 10,
                                         verbose = TRUE)


  imageShow(Region.seg$AP_image_data)



  # change image array to to raster
  sr.1 <- setValues(region.brick, Region.seg$AP_image_data)


  # change to mono colour
  mono.colour <- rgb_2gray(Region.seg$AP_image_data)


  # change to raster
  mon.rast <- raster(mono.colour)
  crs(mon.rast) <- crs(sr.1)
  extent(mon.rast) <- extent(sr.1)

  message("Step 5: k-means clustering starting")

  # kmeans clustering from raster
  vege.df <- raster::as.data.frame(sr.1)

  if(clusters > length(unique(mon.rast))){
    clusters <- length(unique(mon.rast))
  } else   {
    clusters <- clusters
  }

  vege.kmeans <- kmeans(vege.df, clusters)
  k.clusters <- setValues(mon.rast, vege.kmeans$cluster)


  # elbow method
  set.seed(123)
  # function to compute total within-cluster sum of square
  wss <- function(k) {
    kmeans(vege.df, k, nstart = 10 )$tot.withinss
  }

  # Compute and plot wss for k = 1 to max unique clusters
  k.values <- 2:length(unique(mon.rast))

  # extract wss for 2-15 clusters
  wss_values <- map_dbl(k.values, wss)

  elbow.title <- paste0(output.name,"_", type, "_elbow.png")

  elbow.graph <- ggplot()+
    geom_line(aes(x=k.values, y=wss_values))+
    geom_point(aes(x=k.values, y=wss_values))+
    scale_x_continuous(breaks = 1:15)+
    labs(x="\nNumber of clusters",
         y="Total within-clusters sum of squares\n",
         title = paste0(output.name,"_", type, "\nClustering check using elbow method"))

  ggsave(elbow.graph, filename = paste0(out.dir,"\\",  output.name,"_", "_unsupervised_", type, "_elbow.png"))


  message("Step 6: polygon writing")

  # mask sr.1
  final <- mask(k.clusters, region.brick[[1]])

  # change to polygon
  orig.poly <- rasterToPolygons(final, na.rm=TRUE, dissolve = TRUE)

  # change to sf for easier handling
  vege.sf <- st_as_sf(orig.poly)
  vege.sf$ID <- 1:nrow(vege.sf)

  # assign vege.sf to superpixels
  pix.raster <- raster(Region.slic$labels)
  crs(pix.raster) <- crs(sr.1)
  extent(pix.raster) <- extent(sr.1)

  # super pixel rendering
  pix.poly <- rasterToPolygons(pix.raster, na.rm=TRUE, dissolve = TRUE)
  pix.sf <- st_as_sf(pix.poly)



  # write to shp for QGIS
  st_write(vege.sf, paste0(out.dir, "\\", output.name, "_unsupervised_", type,"_clusters", ".shp"),
           driver = "ESRI Shapefile", append =FALSE)
  st_write(pix.sf, paste0(out.dir, "\\", output.name, "_unsupervised_", type,"_superpixel", ".shp"),
           driver = "ESRI Shapefile", append =FALSE)

}


