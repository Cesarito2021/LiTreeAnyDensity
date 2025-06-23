#' Extract XYZ Coordinates and Intensity from Point Cloud Data
#'
#' This function extracts the X, Y, Z coordinates and intensity values from a point cloud data frame and returns them as a matrix.
#' @param data A data frame containing point cloud data with columns named "X", "Y", "Z", and "Intensity".
#' @return A matrix with four columns: X, Y, Z, and Intensity.
#' @examples
#' pc_data <- data.frame(X = runif(10), Y = runif(10), Z = runif(10), Intensity = runif(10))
#' For_xyzi(pc_data)
#' @export
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
#' df <- data.frame(x = runif(100), y = runif(100), z = rnorm(100))
#' r <- ConvertXYZtoRASTER(df, feature = "mean", resolution = 1, crs = "+init=epsg:32632")
#' @export
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
#' @export
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
#' @export
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
#' \dontrun{
#' result <- overlay_func(
#'   data_list,
#'   number_layers = 5,
#'   values_x = rep(1, 5),
#'   values_y = rep(2, 5),
#'   j1 = c(2, 3, 4, 5, 1)
#' )
#' }
#' @export
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

#' Circle Fitting using Sum of Squares Optimization
#'
#' Fits a circle to 2D point data by minimizing the sum of squared differences
#' between observed and predicted radii.
#'
#' @param xy A numeric matrix or data frame with columns representing 2D coordinates (x, y)
#' @param a0 Initial x-coordinate for circle center (default: mean of x values)
#' @param b0 Initial y-coordinate for circle center (default: mean of y values)
#' @param r0 Initial radius estimate (default: mean distance from initial center)
#' @param height Height value associated with the circle (optional)
#'
#' @return A list containing:
#' \describe{
#'   \item{par}{Numeric vector of length 3 containing (a, b, r), where:
#'     \describe{
#'       \item{a}{x-coordinate of fitted circle center}
#'       \item{b}{y-coordinate of fitted circle center}
#'       \item{r}{radius of fitted circle}
#'     }
#'   }
#'   \item{height}{The height value passed to the function (if provided)}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' theta <- runif(50, 0, 2*pi)
#' x <- 5 + 3*cos(theta) + rnorm(50, sd = 0.1)
#' y <- 2 + 3*sin(theta) + rnorm(50, sd = 0.1)
#' xy <- cbind(x, y)
#' result <- fitSS(xy)
#' print(result$par)
fitSS <- function(xy, a0 = mean(xy[, 1]), b0 = mean(xy[, 2]), 
                  r0 = mean(sqrt((xy[, 1] - a0)^2 + (xy[, 2] - b0)^2)), 
                  height = NULL) {
  
  # Input validation
  if (!is.matrix(xy) && !is.data.frame(xy)) {
    stop("xy must be a matrix or data frame")
  }
  if (ncol(xy) < 2) {
    stop("xy must have at least 2 columns for x and y coordinates")
  }
  
  # Sum of squares function to minimize
  SS <- function(abr) {
    sum((abr[3] - sqrt((xy[, 1] - abr[1])^2 + (xy[, 2] - abr[2])^2))^2)
  }
  
  # Optimize using Nelder-Mead method
  res <- optim(c(a0, b0, r0), SS, method = "Nelder-Mead")
  
  # Return results including height if provided
  return(list(par = res$par, height = height))
}

#' Robust Circle Fitting Using the RANSAC Algorithm
#'
#' This function implements the RANSAC (Random Sample Consensus) algorithm to robustly fit a circle to 2D spatial point data (typically x and y coordinates of tree cross-sections), and estimates the average height (z) of the detected shape.
#'
#' @param data A data frame containing at least the columns `x`, `y`, and `z`, representing the 3D coordinates of the points.
#' @param n An integer specifying the minimum number of points required to fit a model (typically 3 for a circle).
#' @param k An integer indicating the maximum number of iterations allowed in the RANSAC algorithm.
#' @param t A numeric threshold value to determine whether a point fits the model (distance tolerance).
#' @param d An integer specifying the minimum number of inliers required to consider a model valid.
#'
#' @return A list representing the best-fit circle model with components `par` (center x, y, and radius), and `height` (mean z-value of the input points).
#' @details The function uses a custom fitting method `fitSS()` to estimate the parameters of a circle from sampled inliers. The model with the lowest average residual error is selected. The mean height (`z`) of all input points is returned as an additional attribute.
#' @export
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

#' Multiple Circle Detection Using the RANSAC Algorithm
#'
#' This function detects and fits multiple circles within a 2D point cloud using the RANSAC (Random Sample Consensus) algorithm. It iteratively applies a robust fitting procedure to identify circular shapes (e.g., tree trunks) and removes inliers after each detection to allow for the next circle to be found.
#'
#' @param data A data frame containing at least the columns `x`, `y`, and `z`, representing the 3D coordinates of the points.
#' @param n An integer specifying the minimum number of points required to fit a model (typically 3 for a circle).
#' @param k An integer indicating the maximum number of iterations allowed in each RANSAC cycle.
#' @param t A numeric threshold that determines whether a point is considered an inlier based on its distance to the circle.
#' @param d An integer specifying the minimum number of inliers required for a circle to be considered valid.
#' @param max_circles Maximum number of circles to detect in the point cloud.
#'
#' @return A list of circle models, each containing the fitted parameters (`par` = center x, center y, radius) and the average height (`height`) of the points that belong to that circle.
#' @details This function calls \code{rc()} internally for each circle to be detected. After each valid detection, inlier points are removed to allow detection of additional non-overlapping circles.
#' @export
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

#' Rotate and Align LAS Point Cloud (Wrapper Function)
#'
#' This function provides a unified interface to rotate and align a 3D mesh point cloud using one of two available methods.
#' The goal is to reorient the point cloud so that the main planar surface (usually the ground or horizontal base) becomes aligned with the vertical Z-axis.
#'
#' @param mesh A mesh object from a `.ply` file, containing a `vb` component with 3D vertex coordinates.
#' @param crs_code Numeric or character. Coordinate Reference System (CRS) code to assign to the output LAS object. Default is EPSG:3035.
#' @param method Integer. Selects the rotation method to apply:
#'   \itemize{
#'     \item \code{1} = PCA-based alignment using the lowest 5% of points (ground-based rotation).
#'     \item \code{2} = PCA alignment of the entire centered point cloud (global variance-based rotation).
#'   }
#'
#' @return A LAS object with rotated and aligned coordinates, and optional CRS assignment.
#'
#' @details
#' \strong{Method 1:} Uses Principal Component Analysis (PCA) on the lowest 5% of Z-values (assumed ground points),
#' then applies Rodrigues' rotation formula to align the normal vector of the estimated plane to the Z-axis.
#'
#' \strong{Method 2:} Centers the full point cloud, performs PCA, and applies the resulting rotation matrix to align the major axis of variation to the Z-axis.
#'
#' Use Method 1 when the base of the object (e.g., ground) is clearly flat and distinguishable.
#' Use Method 2 when the point cloud is more abstract, symmetric, or lacks a clear ground reference.
#'
#' @examples
#' \dontrun{
#' mesh <- vcgPlyRead("path/to/mesh.ply", updateNormals = TRUE)
#' las_rotated <- rotate_and_align_las(mesh, method = 1)
#' }
#'
#' @export
rotate_and_align_las <- function(mesh, crs_code = 3035, method = 1) {
  if (method == 1) {
    rotate_and_align_las_method1(mesh, crs_code)
  } else if (method == 2) {
    rotate_and_align_las_method2(mesh, crs_code)
  } else {
    stop("Invalid method. Please choose method = 1 or method = 2.")
  }
}

#' Rotate and Align LAS Point Cloud (Method 1: Ground-Based PCA)
#'
#' This function rotates a 3D mesh represented by a LAS point cloud object so that its lowest surface points
#' are aligned with the vertical (Z) axis. It applies Principal Component Analysis (PCA) on the lowest 5th percentile
#' of elevation values to estimate the ground normal, then uses Rodrigues' rotation formula to align the point cloud.
#'
#' @param mesh A mesh object containing vertex data, expected to have a component \code{vb} with 3D coordinates.
#' @param crs_code EPSG code for the Coordinate Reference System (CRS) to assign to the LAS object. Default is 3035.
#'
#' @return A LAS object with rotated and vertically aligned coordinates, with CRS and classification fields assigned.
#'
#' @details
#' This method is suited for terrestrial or ground-based point clouds where the lowest points are assumed to represent
#' the ground. It uses PCA on these points to determine the dominant plane and rotates the point cloud accordingly.
#'
#' @examples
#' \dontrun{
#' las_aligned <- rotate_and_align_las_method1(mesh)
#' }
#' @export
rotate_and_align_las_method1 <- function(mesh, crs_code = 3035) {
  rotate_and_align_las_Arch1(mesh, crs_code)
}

#  Note: no @export tag here.
rotate_and_align_las_Arch1 <- function(mesh, crs_code = 3035) {
  vertices <- t(mesh$vb[1:3, ])
  colnames(vertices) <- c("X", "Y", "Z")
  vertices <- as.data.frame(vertices)
  
  suppressMessages(las <- LAS(vertices))
  sf::st_crs(las) <- sf::st_crs(crs_code)
  
  low_points <- las@data[las@data$Z <= quantile(las@data$Z, 0.05), ]
  pca <- prcomp(low_points[, c("X", "Y", "Z")])
  normal_vector <- pca$rotation[, 3]
  
  target <- c(0, 0, 1)
  v <- pracma::cross(normal_vector, target)
  s <- sqrt(sum(v^2))
  c <- sum(normal_vector * target)
  
  vx <- matrix(c(0, -v[3], v[2],
                 v[3], 0, -v[1],
                 -v[2], v[1], 0), 3, 3, byrow = TRUE)
  R <- diag(3) + vx + vx %*% vx * ((1 - c) / (s^2))
  
  coords <- as.matrix(las@data[, c("X", "Y", "Z")])
  coords_rot <- t(R %*% t(coords))
  coords_rot[, 3] <- coords_rot[, 3] - min(coords_rot[, 3])
  
  las@data$X <- coords_rot[, 1]
  las@data$Y <- coords_rot[, 2]
  las@data$Z <- coords_rot[, 3]
  las@data$ReturnNumber <- 1L
  las@data$NumberOfReturns <- 1L
  
  return(las)
}


#' Rotate and Align LAS Point Cloud (Method 2: Global PCA-Based)
#'
#' This function rotates a 3D mesh by applying PCA on the entire centered point cloud to align
#' the axes of variation. It is useful for objects where orientation is not defined by ground points
#' but by overall geometry (e.g., indoor scans or artifacts).
#'
#' @param mesh A mesh object containing vertex data, expected to have a component \code{vb} with 3D coordinates.
#' @param crs_code EPSG code for the Coordinate Reference System (CRS) to assign to the LAS object. Default is 3035.
#'
#' @return A LAS object with globally aligned coordinates using PCA.
#'
#' @details
#' The function reorders axes (X, Z, Y → X, Y, Z), centers the coordinates, performs PCA,
#' and rotates the point cloud using the PCA rotation matrix. This is useful when no clear ground is present.
#'
#' @examples
#' \dontrun{
#' las_aligned <- rotate_and_align_las_method2(mesh)
#' }
#' @export
rotate_and_align_las_method2 <- function(mesh, crs_code = 3035) {
  rotate_and_align_las_Arch2(mesh, crs_code)
}

#  Note: no @export tag here.
rotate_and_align_las_Arch2 <- function(ply, crs_code = 3035) {
  coords <- t(ply$vb[1:3, , drop = FALSE])
  colnames(coords) <- c("X", "Y", "Z")
  coords <- coords[, c("X", "Z", "Y")]
  colnames(coords) <- c("X", "Y", "Z")
  
  coords_centered <- scale(coords, center = TRUE, scale = FALSE)
  points_df <- as.data.table(coords_centered)
  las <- LAS(points_df)
  sf::st_crs(las) <- sf::st_crs(crs_code)
  
  pca <- prcomp(coords_centered)
  rot_mat <- pca$rotation
  rotated_coords <- coords_centered %*% rot_mat
  
  las@data$X <- round(rotated_coords[, 1], 3)
  las@data$Y <- round(rotated_coords[, 2], 3)
  las@data$Z <- round(rotated_coords[, 3], 3)
  
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
#' @export
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
#' @export
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
#' @export
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

#' Slice 3D Point Cloud by Height Intervals
#'
#' This function slices a 3D LiDAR point cloud into multiple horizontal layers within a specified height range and voxel size. It is typically used in forestry applications to analyze structural features such as tree architecture at different vertical strata.
#'
#' @param LasData A LAS object or data frame containing 3D point cloud data, including at least the coordinates `x`, `y`, and `z`.
#' @param Z_mininum A numeric value indicating the minimum height (Z) threshold for filtering the point cloud.
#' @param Z_maximum A numeric value indicating the maximum height (Z) threshold for filtering the point cloud.
#' @param voxel_size A numeric value specifying the resolution (in meters) for both vertical and horizontal slicing (voxel height and width).
#'
#' @return A list with two elements: 
#' \itemize{
#'   \item A list of data frames representing the sliced horizontal layers of the point cloud.
#'   \item A data frame containing the filtered `x`, `y`, and `z` coordinates used for slicing.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming `las` is a loaded LAS file object
#' result <- Slice3DByHeightArch(las, Z_mininum = 2, Z_maximum = 20, voxel_size = 0.5)
#' slices <- result[[1]]
#' xyz_data <- result[[2]]
#' }
#  Note: no @export tag here.
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

#' Rasterization of 3D Sliced Point Cloud Data
#'
#' This function rasterizes a list of 3D point cloud slices (e.g., horizontal layers of a tree canopy)
#' and computes spatial metrics such as the number of occupied voxels and the sum of values across slices.
#' It is designed for use in forest structure analysis from LiDAR data.
#'
#' @param LasData A list of data frames, each containing 3D point cloud data (columns: `x`, `y`, `z`) 
#' representing slices of the original point cloud (e.g., from \code{Slice3DByHeightArch}).
#' @param voxel_size A numeric value specifying the spatial resolution (in meters) of the output rasters.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{RasterMetric}: A list of binary rasters indicating occupied voxels in each slice.
#'   \item \code{RasterMetricSum}: A raster showing the total number of layers each pixel was occupied.
#'   \item \code{RasterMetricCount}: A raster showing the sum of point values across all slices.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'sliced_data' is the output from Slice3DByHeightArch
#' voxel_results <- Rasterize3DSlicesArch(sliced_data[[1]], voxel_size = 0.5)
#' plot(voxel_results[[2]])  # RasterMetricSum
#' }
#  Note: no @export tag here.
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

#' Identification of Potential Tree Stem Locations from Rasterized Point Cloud Metrics
#'
#' This function identifies potential tree stem positions from a raster stack derived from 3D point cloud slices.
#' It evaluates voxel metrics across multiple canopy layers to detect persistent vertical structures, 
#' applying spatial filtering and thresholding to locate likely stem base positions.
#'
#' @param xyz_data A data frame with columns \code{x}, \code{y}, and \code{z}, representing 3D coordinates of the original point cloud.
#' @param RasterMetric A list of binary rasters indicating voxel occupancy for each canopy layer (e.g., from \code{Rasterize3DSlicesArch}).
#' @param RasterMetricSum A raster representing the number of occupied layers per pixel.
#' @param RasterMetricCount A raster representing the total point count across all layers per pixel.
#' @param voxel_size A numeric value indicating the voxel resolution used to compute the raster metrics (in meters).
#'
#' @return A data frame of filtered 3D points representing potential tree stem base locations.
#'
#' @examples
#' \dontrun{
#'result <- IdentifyPotentialTreeLocationsArch(
#'  xyz_data, RasterMetric, RasterMetricSum, RasterMetricCount,
#'  voxel_size = 0.5)
#'}
#'
#  Note: no @export tag here.
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

#' Tree Detection and Diameter Measurement from Point Cloud Data
#'
#' This function detects tree stems and estimates their diameters using DBSCAN clustering followed by RANSAC-based circle fitting on 3D point cloud data. It outputs an `sf` object with the estimated diameter at breast height (DBH), basal area, and stem location. This tool is designed for use in forest structure analysis from terrestrial LiDAR or photogrammetric point clouds.
#'
#' @param data_ref A data frame containing 3D coordinates (x, y, z) of point cloud data.
#' @param eps_value The epsilon parameter for the DBSCAN clustering algorithm, defining the neighborhood radius.
#' @param minPts_value The minimum number of points required to form a dense region in DBSCAN.
#' @param d1_m_minimum Minimum accepted DBH (in meters) for filtering detected stems.
#' @param d1_m_maximum Maximum accepted DBH (in meters) for filtering detected stems.
#' @param n_ransac_par The minimum number of points in a cluster required to apply RANSAC.
#' @param k_ransac_par The number of RANSAC iterations.
#' @param t_ransac_par The threshold distance to determine inliers during RANSAC fitting.
#' @param d_ransac_par The distance threshold for merging or validating detected shapes.
#'
#' @return An `sf` object (EPSG:3035) containing polygons of detected stem cross-sections with columns: `id` (cluster ID), `dbh_cm` (diameter in cm), `g_m2` (basal area in m²), and `geometry` (buffered stem cross-section).
#'
#' @examples
#' \dontrun{
#' data_ref <- data.frame(x = runif(100), y = runif(100), z = runif(100))
#' result <- DetectAndMeasureTreesArch(data_ref, eps_value = 0.3, minPts_value = 10,
#'                                     d1_m_minimum = 0.05, d1_m_maximum = 1.0,
#'                                     n_ransac_par = 30, k_ransac_par = 100,
#'                                     t_ransac_par = 0.01, d_ransac_par = 0.05)
#' plot(result["dbh_cm"])
#' }
#'
#  Note: no @export tag here.
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


