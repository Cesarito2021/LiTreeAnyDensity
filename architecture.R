#' Extract XYZ Coordinates and Intensity from Point Cloud Data
#'
#' This function extracts the X, Y, Z coordinates and intensity values from a point cloud data frame and returns them as a matrix.
#' @param data A data frame containing point cloud data with columns named "X", "Y", "Z", and "Intensity".
#' @return A matrix with four columns: X, Y, Z, and Intensity.
#' @examples
#' pc_data <- data.frame(X = runif(10), Y = runif(10), Z = runif(10), Intensity = runif(10))
#' For_xyzi(pc_data)
For_xyzi <- function(data){
  x <- data[["X"]]
  y <- data[["Y"]]
  z <- data[["Z"]]
  i <- data[["Intensity"]]
  return(cbind(x,y,z,i))
}
#' Convert Point Cloud XYZ Data to Raster
#'
#' This function converts point cloud data (x, y, z) into a raster layer based on a specified feature, such as the mean elevation or point count per cell.
#' @param data A data frame containing point cloud data with columns named "x", "y", and "z".
#' @param feature A character string specifying the feature to compute. Options are \code{"mean"} for mean z-values or \code{"count"} for point density.
#' @param resolution A numeric value specifying the resolution (pixel size) of the output raster.
#' @param crs A valid coordinate reference system (CRS) object or string (e.g., EPSG code) to assign to the output raster.
#' @return A \code{RasterLayer} object representing the spatial distribution of the selected feature.
#' @examples
#' library(raster)
#' df <- data.frame(x = runif(100), y = runif(100), z = rnorm(100))
#' r <- ConvertXYZtoRASTER(df, feature = "mean", resolution = 1, crs = "+init=epsg:32632")
ConvertXYZtoRASTER <- function(data, feature, resolution, crs) {
  if (feature == "mean") {
    ext <- raster::extent(data[, c("x", "y")])
    ext <- ext + as.numeric(resolution)  # Aggiunta del buffer se necessario
    pixel_size <- as.numeric(resolution)  # Specifica la dimensione desiderata dei pixel nelle unità del sistema di coordinate
    r <- raster::raster(ext, resolution = pixel_size)
    means <- raster::rasterize(raster::as.data.frame(data[, c("x", "y")]), r, data[, "z"], fun = mean)
    raster::crs(means) <- crs
    return(means)
  } else if (feature == "count") {
    ext <- raster::extent(data[, c("x", "y")])
    ext <- ext + as.numeric(resolution)  # Aggiunta del buffer se necessario
    pixel_size <- as.numeric(resolution)  # Specifica la dimensione desiderata dei pixel nelle unità del sistema di coordinate
    r <- raster::raster(ext, resolution = pixel_size)
    count <- raster::rasterize(raster::as.data.frame(data[, c("x", "y")]), r, fun = "count")
    raster::crs(count) <- crs
    return(count)
  }
}
# step 2
#' Align Raster Extent and Resolution to a Reference Raster
#'
#' This function adjusts the extent, resolution, and coordinate reference system (CRS) of a raster object to match a reference raster.
#' @param raster_ref A \code{RasterLayer} object used as the reference for extent, resolution, and CRS.
#' @param raster_to_correct A \code{RasterLayer} object that will be resampled and aligned to match the reference raster.
#' @return A \code{RasterLayer} object with extent, resolution, and CRS aligned to \code{raster_ref}.
#' @examples
#' library(raster)
#' r1 <- raster(ncol=10, nrow=10)
#' extent(r1) <- c(0, 10, 0, 10)
#' crs(r1) <- CRS("+init=epsg:4326")
#' r2 <- raster(ncol=5, nrow=5)
#' extent(r2) <- c(1, 6, 1, 6)
#' crs(r2) <- CRS("+init=epsg:3857")
#' r_corrected <- Extent_validation(r1, r2)
Extent_validation <- function(raster_ref,raster_to_correct){
  #
  ref_extent <- extent(raster_ref)
  # 
  raster_to_correct <- resample(raster_to_correct, raster(ref_extent, ncols = ncol(raster_to_correct), nrows = nrow(raster_to_correct)))
  # 
  crs(raster_to_correct) <- crs(raster_ref)
  #
  raster_corrected  <- resample(raster_to_correct, raster_ref, method = "bilinear")
  #
  return(raster_corrected)
}
#' Custom Vectorized Overlay Function
#'
#' This function performs a conditional check on two raster layer values and returns a modified value based on the comparison with given constants.
#'
#' @param x A numeric vector representing values from the first raster layer.
#' @param y A numeric vector representing values from the second raster layer.
#' @param valx A numeric value to compare against elements in \code{x}.
#' @param valy A numeric value to compare against elements in \code{y}.
#'
#' @return A numeric vector where elements are \code{x + 1} if \code{x == valx} and \code{y == valy}, otherwise 0.
#'
#' @details This function is designed to be used within raster overlay operations. It is vectorized and can be applied efficiently to raster data using functions like \code{\link[raster]{overlay}}.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' y <- c(5, 2, 3, 4)
#' vectorized_custom_fun_bis(x, y, valx = 3, valy = 3)
#' # Returns: 0 0 4 0
vectorized_custom_fun_bis <- function(x, y, valx, valy) {
  return(ifelse(x == as.numeric(valx) & y == as.numeric(valy), x + 1, 0))
}
#' Overlay Multiple Raster Layers with Custom Function
#'
#' This function performs pairwise overlay operations on a list of raster layers using a custom vectorized function. Each overlay operation applies specific values from \code{values_x} and \code{values_y} to modify the overlay result.
#'
#' @param data A list of \code{RasterLayer} objects to be processed.
#' @param number_layers An integer specifying the number of layers to process.
#' @param values_x A numeric vector of values to be passed to the custom overlay function.
#' @param values_y A numeric vector of values to be passed to the custom overlay function.
#' @param j1 An integer vector of indices indicating which layers in \code{data} should be overlaid with each corresponding layer.
#'
#' @return A list of \code{RasterLayer} objects, each resulting from the custom overlay operation.
#'
#' @details The function assumes the existence of an external function named \code{vectorized_custom_fun_bis}, which defines the overlay behavior between two raster layers and uses the parameters \code{valx} and \code{valy}.
#'
#' @examples
#' # Example not run because it requires pre-defined raster layers and vectorized_custom_fun_bis
#' # result <- overlay_func(data_list, number_layers = 5, values_x = rep(1, 5), values_y = rep(2, 5), j1 = c(2, 3, 4, 5, 1))
overlay_func <- function(data, number_layers, values_x, values_y, j1) {
  output_list <- vector("list", number_layers)
  
  for (i in 1:min(length(data), length(values_x))) {
    valx <- values_x[i]
    valy <- values_y[i]
    j2 <- j1[i]
    
    if (!is.null(data[[i]]) && !is.null(data[[j2]])) {
      if (is.null(output_list[[i]])) {
        output_list[[i]] <- overlay(data[[i]], data[[j2]], fun = function(x, y) vectorized_custom_fun_bis(x, y, valx, valy))
      } else {
        output_list[[i]] <- overlay(output_list[[i]], data[[j2]], fun = function(x, y) vectorized_custom_fun_bis(x, y, valx, valy))
      }
    } else {
      print(paste("Skipped due to NULL data at index:", i))
    }
  }
  
  return(output_list)
}
#' Circle Fitting and Multiple Circle Detection in 2D Point Data
#'
#' These functions implement circle fitting and robust circle detection using a RANSAC-like approach.
#'
#' @param xy A numeric matrix or data frame with columns representing 2D coordinates (x, y) for circle fitting (used in \code{fitSS}).
#' @param data A data frame or matrix with columns \code{x}, \code{y}, and \code{z} representing 3D points (used in \code{rc} and \code{rmc}).
#' @param n Integer, minimum number of points used to fit a circle in each iteration.
#' @param k Integer, maximum number of iterations for circle detection.
#' @param t Numeric, threshold distance to consider a point as an inlier to a circle.
#' @param d Integer, minimum number of inliers required to accept a circle model.
#' @param max_circles Integer, maximum number of circles to detect (used in \code{rmc}).
#'
#' @return
#' \code{fitSS} returns a list with estimated circle parameters (center coordinates and radius).
#' \code{rc} returns the best fitted circle model found after RANSAC iterations, including the average height.
#' \code{rmc} returns a list of fitted circle models detected in the data.
#'
#' @details
#' \code{fitSS} fits a single circle to 2D points by minimizing squared residuals.
#' \code{rc} applies a RANSAC-like method to find the best circle fitting a subset of points.
#' \code{rmc} iteratively detects multiple circles by applying \code{rc} and removing inliers.
#'
#' @examples
#' \dontrun{
#'   single_circle <- fitSS(matrix(c(x, y), ncol=2))
#'   best_circle <- rc(data, n=10, k=100, t=0.1, d=15)
#'   circles <- rmc(data, n=10, k=100, t=0.1, d=15, max_circles=5)
#' }
fitSS <- function(xy, a0 = mean(xy[, 1]), b0 = mean(xy[, 2]), r0 = mean(sqrt((xy[, 1] - a0)^2 + (xy[, 2] - b0)^2)), height) {
  SS <- function(abr) {
    sum((abr[3] - sqrt((xy[, 1] - abr[1])^2 + (xy[, 2] - abr[2])^2))^2)
  }
  
  res <- optim(c(a0, b0, r0), SS, method = "Nelder-Mead")  # Cambiato L-BFGS-B a Nelder-Mead
  
  return(list(par = res$par))
}
rc <- function(data, n, k, t, d) {
  bestfit <- NULL
  besterr <- Inf
  
  for (iterations in 1:k) {
    maybeinliers <- sample(nrow(data), n)
    maybemodel <- fitSS(data[maybeinliers, c("x", "y")])
    
    distances <- sqrt((data[, "x"] - maybemodel$par[1])^2 + (data[, "y"] - maybemodel$par[2])^2)
    inliers <- which(abs(distances - maybemodel$par[3]) < t)
    
    if (length(inliers) > d) {
      bettermodel <- fitSS(data[c(maybeinliers, inliers), c("x", "y")])
      residuals <- abs(distances[inliers] - bettermodel$par[3])
      thiserr <- mean(residuals)  # Cambiato sd a mean per la velocità
      
      if (thiserr < besterr) {
        bestfit <- bettermodel
        besterr <- thiserr
      }
    }
  }
  
  bestfit$height <- mean(data[, "z"])  # Altezza media
  return(bestfit)
}
rmc <- function(data, n, k, t, d, max_circles) {
  all_circles <- vector("list", max_circles)
  remaining_points <- data
  found_circles <- 0
  
  while (nrow(remaining_points) >= n && found_circles < max_circles) {
    bestfit <- rc(remaining_points, n, k, t, d)
    
    if (!is.null(bestfit)) {
      all_circles[[found_circles + 1]] <- bestfit
      
      distances <- sqrt((remaining_points[, "x"] - bestfit$par[1])^2 + (remaining_points[, "y"] - bestfit$par[2])^2)
      inliers <- which(abs(distances - bestfit$par[3]) < t)
      
      remaining_points <- remaining_points[-inliers, , drop = FALSE]  # Rimuovi gli inliers
      found_circles <- found_circles + 1
    } else {
      break
    }
  }
  
  return(all_circles[1:found_circles])
}
#
#' Rotate and Align LAS Point Cloud to Vertical Axis
#'
#' This function rotates a 3D mesh represented by a LAS point cloud object so that its lowest surface points
#' are aligned with the vertical (Z) axis. It applies Principal Component Analysis (PCA) on the lowest 5th percentile
#' of points to estimate the normal vector, then uses Rodrigues' rotation formula to realign the mesh.
#'
#' @param mesh A mesh object containing vertex data, expected to have a component \code{vb} with 3D coordinates.
#' @param crs_code Numeric or character representing the Coordinate Reference System (CRS) code to assign to the LAS object. Default is EPSG:3035.
#'
#' @return A LAS object with rotated and vertically aligned coordinates, with updated CRS and basic classification fields.
#'
#' @details
#' The function extracts vertex coordinates from the mesh, constructs a LAS object, estimates the ground plane via PCA on
#' the lowest 5% Z values, calculates the rotation matrix to align the ground plane normal to the vertical axis,
#' applies the rotation, and repositions the point cloud so the minimum Z is zero.
#' 
#' @examples
#' \dontrun{
#'   las_aligned <- rotate_and_align_las(mesh_object)
#' }
rotate_and_align_las <- function(mesh, crs_code = 3035) {
  # 1. Estrai vertici dal mesh e crea data.frame
  vertices <- t(mesh$vb[1:3, ])
  colnames(vertices) <- c("X", "Y", "Z")
  vertices <- as.data.frame(vertices)
  
  # 2. Crea LAS e assegna CRS
  suppressMessages(las <- LAS(vertices))
  sf::st_crs(las) <- sf::st_crs(crs_code)
  
  # 3. Seleziona punti bassi (5° percentile Z) per PCA
  low_points <- las@data[las@data$Z <= quantile(las@data$Z, 0.05), ]
  pca <- prcomp(low_points[, c("X","Y","Z")])
  normal_vector <- pca$rotation[,3]
  
  # 4. Calcola matrice di rotazione (Rodrigues) per allineare normale a vettore Z
  target <- c(0,0,1)
  v <- pracma::cross(normal_vector, target)
  s <- sqrt(sum(v^2))
  c <- sum(normal_vector * target)
  vx <- matrix(c(0, -v[3], v[2],
                 v[3], 0, -v[1],
                 -v[2], v[1], 0), 3, 3, byrow = TRUE)
  R <- diag(3) + vx + vx %*% vx * ((1 - c)/(s^2))
  
  # 5. Applica rotazione e riallinea Z a 0
  coords <- as.matrix(las@data[, c("X","Y","Z")])
  coords_rot <- t(R %*% t(coords))
  coords_rot[,3] <- coords_rot[,3] - min(coords_rot[,3])
  
  # 6. Aggiorna dati LAS
  las@data$X <- coords_rot[,1]
  las@data$Y <- coords_rot[,2]
  las@data$Z <- coords_rot[,3]
  
  # 7. Imposta classificazione base
  las@data$ReturnNumber <- 1L
  las@data$NumberOfReturns <- 1L
  
  return(las)
}
#' Filter LAS Points Within the Central Area
#'
#' This function filters points in a LAS object to retain only those located within a central area,
#' excluding a specified margin percentage from the edges in the XY plane.
#'
#' @param las A LAS object containing 3D point cloud data with X and Y coordinates.
#' @param margin_pct A numeric value between 0 and 0.5 indicating the percentage of the range to exclude from each side (default is 0.05).
#'
#' @return A LAS object containing only the points within the central XY area defined by the margin.
#'
#' @details
#' The function calculates the XY extent of the LAS points, defines a reduced central area by removing margins
#' proportional to \code{margin_pct} from all sides, and returns the subset of points falling within this area.
#'
#' @examples
#' \dontrun{
#'   las_center <- filter_las_central_area(las_object, margin_pct = 0.1)
#' }
filter_las_central_area <- function(las, margin_pct = 0.05) {
  # Calcola limiti XY
  xrange <- range(las@data$X)
  yrange <- range(las@data$Y)
  
  # Calcola margini
  x_margin <- margin_pct * diff(xrange)
  y_margin <- margin_pct * diff(yrange)
  
  # Definisci area centrale
  xmin <- xrange[1] + x_margin
  xmax <- xrange[2] - x_margin
  ymin <- yrange[1] + y_margin
  ymax <- yrange[2] - y_margin
  
  # Filtra punti nell'area centrale
  las_filtered_center <- filter_poi(las, X >= xmin & X <= xmax & Y >= ymin & Y <= ymax)
  
  return(las_filtered_center)
}
#' Create Circle Attributes from RANSAC Output
#'
#' This function generates a data frame with geometric and identifier attributes for a detected circle
#' based on parameters obtained from a RANSAC fitting output. It applies optional diameter filtering
#' and prepares the data for spatial processing.
#'
#' @param ransac_output A list containing circle parameters (`par`), typically from a RANSAC circle fitting function.
#' @param circle_id An integer or character identifier for the circle.
#' @param crs Coordinate Reference System (CRS) code (EPSG) to assign to the output (currently not used).
#' @param d1_m_minimum Minimum allowed diameter for the circle (numeric). Circles with diameter below this are discarded.
#' @param d1_m_maximum Maximum allowed diameter for the circle (numeric). Circles with diameter above this are discarded.
#' @param n_points Number of points supporting the circle (currently not used).
#'
#' @return A data frame with the circle model type, id, centroid coordinates (`centroid_x`, `centroid_y`),
#' diameter (`d1_m`), and area (`area_m2`). Returns NULL if parameters are missing or diameter filters fail.
#'
#' @examples
#' \dontrun{
#'  result <- create_circle_shapefile_from_ransac(ransac_output = list(par = c(10, 20, 5)),
#'                                               circle_id = 1,
#'                                               crs = 32633,
#'                                               d1_m_minimum = 2,
#'                                               d1_m_maximum = 10,
#'                                               n_points = 50)
#' }
create_circle_shapefile_from_ransac <- function(ransac_output, circle_id, crs, d1_m_minimum, d1_m_maximum, n_points) {
  # Check if the required elements are present
  if (is.null(ransac_output$par) || length(ransac_output$par) < 3) {
    warning(paste("RANSAC output for circle", circle_id, "is missing parameters."))
    return(NULL)
  }
  # Extract center coordinates and radius from `par`
  center_x <- ransac_output$par[1]  # First parameter
  center_y <- ransac_output$par[2]  # Second parameter
  radius <- ransac_output$par[3]     # Third parameter (assuming it represents the radius)
  dbh <- round(radius, 2)
  ba <-pi*(radius)^2
  id <- circle_id
  #
  # Apply diameter filters
  if (!is.null(d1_m_minimum) && !is.null(d1_m_maximum)) {
    if (dbh < d1_m_minimum  || dbh > d1_m_maximum ) {
      warning("Circle diameter does not meet the thresholds.")
      return(NULL)
    }
  }
  # Create a properties table
  tab_circle <- data.frame(
    Model = "Circle",
    id = id,
    centroid_x = center_x,
    centroid_y = center_y,
    d1_m = round(as.numeric(dbh), 2),  # Round to 2 decimals
    area_m2 = round(as.numeric(ba), 2)  # Round to 2 decimals
  )
  #
  #tab_circle_sf <- st_as_sf(tab_circle,coords = c('centroid_x','centroid_y'),crs=32633)
  #tab_circle_sf$geometry <- st_buffer(tab_circle_sf$geometry, dist=tab_circle_sf$d1_m/2)
  return(tab_circle)
}
#' Calculate Cluster Centroids and Assign Cluster IDs
#'
#' This function computes the centroid coordinates (mean of x, y, and z) for each cluster in a dataset
#' and assigns cluster IDs to the individual points. It returns a data frame containing the original points,
#' their assigned cluster IDs, and the corresponding cluster centroids.
#'
#' @param data A data frame or matrix with point coordinates, typically including columns for x, y, and z.
#' @param cluster_labels A vector of cluster labels (integers or factors) corresponding to each row in `data`.
#'
#' @return A data frame including the original point coordinates, a column `ID_Cluster` with cluster assignments,
#' and two additional columns `Centroid_X` and `Centroid_Y` representing the centroid coordinates of each cluster.
#'
#' @examples
#' \dontrun{
#'  pts <- data.frame(x = runif(10), y = runif(10), z = runif(10))
#'  labels <- c(1,1,2,2,1,3,3,2,1,3)
#'  result <- calculate_cluster_centroids(pts, labels)
#' }

calculate_cluster_centroids <- function(data, cluster_labels) {
  # Calculate centroids
  centroids <- lapply(unique(cluster_labels), function(label) {
    cluster_points <- data[cluster_labels == label, ]
    colMeans(cluster_points)
  })
  
  # Create an empty table for cluster points
  cluster_points_table <- data.frame(x = numeric(),
                                     y = numeric(),
                                     z = numeric(),
                                     ID_Cluster = integer())
  
  # Iterate through the clusters
  for (label in unique(cluster_labels)) {
    # Get points of current cluster
    cluster_points <- data[cluster_labels == label, ]
    
    # Calculate centroid of current cluster
    centroid <- centroids[[label]]
    
    # Add cluster ID to each row
    cluster_points <- dplyr::mutate(cluster_points, ID_Cluster = label)
    
    # Add current cluster points to the table
    cluster_points_table <- bind_rows(cluster_points_table, cluster_points)
  }
  
  # Merge centroids with the cluster points table
  cluster_points_table <- left_join(cluster_points_table, data.frame(ID_Cluster = unique(cluster_labels), 
                                                                     Centroid_X = sapply(centroids, "[[", 1),
                                                                     Centroid_Y = sapply(centroids, "[[", 2)), by = "ID_Cluster")
  
  return(cluster_points_table)
}


Slice3DByHeightArch <- function(LasData,Z_mininum,Z_maximum,voxel_size){
  # cut the horizontal slices
  filtered_tree_slice<- filter_poi(LasData, Z>Z_mininum & Z <Z_maximum)
  # convert las file in xyz
  data_ref <- as.data.frame(For_xyzi(filtered_tree_slice))
  # extract xyz data
  xyz_data<-  data_ref %>% dplyr::select(c("x","y","z"))
  #--- set voxel height ---#
  vs=voxel_size
  #--- set voxel horizontal ---#
  vh=voxel_size
  # set number of voxel layers for finding trees
  vn=10
  # voxel analysis in the height level between "lower_limit" and "upper_limit" as basis for detecting trees
  lower_limit = 1
  upper_limit = lower_limit+vh*vn 
  #-------------------------------------------------
  #-------------- Extract A Slice Horizintal -------
  #-------------------------------------------------
  number_slices = 10
  #
  multi_slices  <-  xyz_data %>%dplyr::mutate(interval = cut(z, breaks = seq(min(z), max(z),c(max(z)-min(z))/number_slices), include.lowest = TRUE))
  #
  LasData  <- multi_slices%>% group_split(interval)
  #
  return(list(LasData,xyz_data))}
Rasterize3DSlicesArch <- function(LasData,voxel_size){
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------
  m <- foreach(i = 1:length(LasData)) %do% {
    ConvertXYZtoRASTER(LasData[[i]], feature = "mean", resolution = voxel_size, crs = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
  }
  #-------------------------------------------------
  m1 <- suppressWarnings(
    foreach(i = 1:length(m)) %do% {
      Extent_validation(raster_ref = m[[1]],raster_to_correct = m[[i]])})
  # Crea lo stack di raster utilizzando la funzione stack
  stacked_raster <- do.call(stack, m1)
  # Calcola la somma dei raster utilizzando la funzione 'calc'
  RasterMetricCount<- raster::calc(stacked_raster, sum,na.rm=T)
  #
  RasterMetric  <- foreach(i = 1:length(m1)) %do% {
    raster::calc(m1[[i]], function(x) ifelse(x > 0, 1, NA))}  
  # Crea lo stack di raster utilizzando la funzione stack
  stacked_raster_p1b_p10b <- do.call(stack, RasterMetric)
  # Calcola la somma dei raster utilizzando la funzione 'calc'
  RasterMetricSum <-  calc(stacked_raster_p1b_p10b, sum,na.rm=T)
  #
  return(list(RasterMetric,RasterMetricSum,RasterMetricCount))
}
IdentifyPotentialTreeLocationsArch  <- function(xyz_data, RasterMetric, RasterMetricSum, RasterMetricCount, voxel_size) {
  # Define values for each case
  layers_params <- list(
    list(number_layers = 10, values_x = 1:10, values_y = rep(1, 10), j1 = c(2, 3, 4, 5, 6, 7, 8, 9, 9, 10)),
    list(number_layers = 9, values_x = 1:9, values_y = rep(1, 9), j1 = c(2, 3, 4, 5, 6, 7, 8, 9, 9)),
    list(number_layers = 8, values_x = 1:8, values_y = rep(1, 8), j1 = c(2, 3, 4, 5, 6, 7, 8, 9)),
    list(number_layers = 7, values_x = 1:7, values_y = rep(1, 7), j1 = c(2, 3, 4, 5, 6, 7, 8)),
    list(number_layers = 6, values_x = 1:6, values_y = rep(1, 6), j1 = c(2, 3, 4, 5, 6, 7)),
    list(number_layers = 5, values_x = 1:5, values_y = rep(1, 5), j1 = c(2, 3, 4, 5, 6)),
    list(number_layers = 4, values_x = 1:4, values_y = rep(1, 4), j1 = c(2, 3, 4, 5)),
    list(number_layers = 3, values_x = 1:3, values_y = rep(1, 3), j1 = c(2, 3, 4)),
    list(number_layers = 2, values_x = 1:2, values_y = rep(1, 2), j1 = c(2, 3))
  )
  
  # Apply overlay function for each set of parameters
  output_list_all <- lapply(layers_params, function(params) {
    print(paste("Processing layers with", params$number_layers, "layers"))
    overlay_func(data = RasterMetric, 
                 number_layers = params$number_layers, 
                 values_x = params$values_x, 
                 values_y = params$values_y, 
                 j1 = params$j1)
  })
  
  # Stack and process the layers
  stacked_rasters <- lapply(output_list_all, function(output_list) {
    do.call(stack, output_list)
  })
  
  # Calculate the max layer for each set
  max_layers <- lapply(stacked_rasters, function(stacked_raster) {
    calc(stacked_raster, max)
  })
  
  # Combine all the max layers if needed (e.g., for the final output)
  final_stack <- do.call(stack, max_layers)
  
  # Remove files with this pattern
  kernelSize <- 1
  max_layers_count <- calc(final_stack, max)
  max_layers_count_st3 <- focal(max_layers_count, w = matrix(1, nrow = (2 * kernelSize + 1), ncol = (2 * kernelSize + 1)), fun = mean, na.rm = TRUE)
  
  # Custom function to find potential stems
  ToFindPotentialStems <- function(a, b, c, d) {
    ifelse(a > 1 & b > 0.051 & c > 5, 1, 0)
  }
  
  # Apply overlay function
  v1_25_pot_stem <- overlay(max_layers_count, max_layers_count_st3, RasterMetricCount, RasterMetricSum, fun = ToFindPotentialStems)
  
  # Add stem position to XYZ.txt
  values <- raster::extract(v1_25_pot_stem, cbind(xyz_data["x"], xyz_data["y"]))
  xyz_data_v2 <- xyz_data %>% dplyr::mutate(stem_cl = values)
  #
  # Filter data
  data_ref <- xyz_data_v2 %>% dplyr::filter(stem_cl == 1) %>% dplyr::select(x:z)
  
  # Return filtered data
  return(data_ref)
}
DetectAndMeasureTreesArch <- function(data_ref,eps_value,minPts_value,d1_m_minimum,
                                      d1_m_maximum, n_ransac_par, k_ransac_par,
                                      t_ransac_par, d_ransac_par){
  
  #
  if (is.null(data_ref) || nrow(data_ref) == 0) {
    cat("Empty or NULL list for data_ref\n")
    return(NULL)
  }
  #
  DBCAN_max <-suppressWarnings( dbscan(data_ref, eps=eps_value, minPts = minPts_value))
  # 
  cluster_labels  <- DBCAN_max[["cluster"]]
  #
  cluster_labels_v2   <-ifelse( cluster_labels %in% 0,max(cluster_labels)+1,cluster_labels)
  #
  cluster_points_table <- calculate_cluster_centroids(data = data_ref, cluster_labels = cluster_labels_v2)
  # 
  xyid <- cluster_points_table %>%
    dplyr::select(c(x, y, z, ID_Cluster)) %>% 
    rename(x = x, y = y, id = ID_Cluster)
  #
  xyid_filtrato <- xyid%>%
    group_by(id) %>%
    filter(n() > n_ransac_par) %>%
    ungroup()
  #
  xyidl <- xyid_filtrato %>%group_split(id) 
  # 
  
  xyidl_df <- suppressMessages (
    foreach(j = seq_along(xyidl)) %do% {
      print(paste0("Convert data in data.frame: ", j))
      xyidl[[j]]%>%data.frame()
    }
  )
  #
  results <- suppressMessages (
    foreach(j = seq_along(xyidl_df)) %do% {
      print(paste0("Detecting the diameter of the cluster: ", j))
      rmc(xyidl_df[[j]],n = n_ransac_par, k = k_ransac_par, t = t_ransac_par, d = d_ransac_par, max_circles = 1)  # max_cylinders = 3 è un esempio
    }
  )
  #
  results1      <- Filter(Negate(is.null),results)
  #
  flattened_results <- unlist(lapply(results1, function(x) {
    if (length(x) > 0) {
      return(x)  
    }
  }), recursive = FALSE)
  #
  suppressWarnings(results_dbh <- foreach(k = seq_along(flattened_results)) %do% {
    print(paste0("Measuring-Refining the diameter of the cluster: ", k))
    create_circle_shapefile_from_ransac(flattened_results[[k]], 
                                        circle_id = k, 
                                        crs = 3035, 
                                        d1_m_minimum = d1_m_minimum, 
                                        d1_m_maximum = d1_m_maximum, 
                                        n_points=180)  # Aggiungi n_points
  })
  # # filter negative files
  results_dbh <- Filter(Negate(is.null), results_dbh)
  # # merged of files
  merged_shapefile <- bind_rows(results_dbh)
  # # Basal area calculation
  merged_shapefile <- st_as_sf(merged_shapefile, coords = c("centroid_x", "centroid_y"), crs = 3035)
  merged_shapefile <- merged_shapefile %>%
    mutate(
      dbh_cm = round(d1_m * 100,2),
      g_m2 = round((pi * (d1_m )^2)/4,2)
    )
  # Buffer individuale per ogni geometria
  library(purrr)
  crs_val <- as.numeric(3035)
  merged_shapefile$geometry <- purrr::map2(
    merged_shapefile$geometry,
    merged_shapefile$d1_m / 2,
    ~ st_buffer(.x, dist = .y)
  ) %>% st_sfc(crs = crs_val)
  #
  # Seleziona e rinomina le colonne desiderate
  merged_shapefile <- merged_shapefile |>
    dplyr::select(
      id = id,     # o il nome della colonna che rappresenta l'ID
      dbh_cm,
      g_m2,
      geometry
    )
  
  return(merged_shapefile)
}
