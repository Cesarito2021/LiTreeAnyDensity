# *********************************
# General Functions
# *********************************
#' Extract 3D Horizontal Slices from a LAS Point Cloud
#'
#' This function filters a LAS point cloud within a given height range (`Z_mininum` to `Z_maximum`) and 
#' divides it into evenly spaced horizontal slices. It returns a list containing the point cloud 
#' split into slices and the original (filtered) XYZ coordinates. It is designed to facilitate 
#' voxel-based analysis or tree detection within defined height intervals.
#'
#' @param LasData A LAS object containing 3D point cloud data.
#' @param Z_mininum Numeric. The minimum height threshold (in meters) for slicing the point cloud.
#' @param Z_maximum Numeric. The maximum height threshold (in meters) for slicing the point cloud.
#' @param voxel_size Numeric. The voxel size (in meters) used to define the vertical resolution of slices.
#'
#' @return A list of two elements: 
#' \describe{
#'   \item{\code{LasData}}{A list of data frames, each corresponding to a height slice of the point cloud.}
#'   \item{\code{xyz_data}}{A data frame of filtered XYZ coordinates within the specified height range.}
#' }
#'
#' @examples
#' \dontrun{
#'   library(lidR)
#'   las <- readLAS("example.las")
#'   slices <- Slice3DByHeight(las, Z_mininum = 1, Z_maximum = 10, voxel_size = 0.2)
#'   print(length(slices[[1]])) # number of horizontal slices
#'   head(slices[[2]])          # XYZ coordinates in filtered height range
#' }
#'
#' @export
Slice3DByHeight <- function(LasData,Z_mininum,Z_maximum,voxel_size) {
  Slice3DByHeightArch(LasData,Z_mininum,Z_maximum,voxel_size)
}

#' Rasterization of 3D Point Cloud Slices
#'
#' This function converts a list of 3D point cloud slices into raster layers using a specified voxel size,
#' aligns their extents, and computes raster-based metrics. It is typically used in forest structure analysis,
#' voxelization workflows, or for converting 3D LiDAR data into 2D raster layers for further analysis.
#'
#' @param LasData A list of data frames, each containing XYZ coordinates from a 3D point cloud slice. 
#' Typically the output from \code{\link{Slice3DByHeight}}.
#' @param voxel_size Numeric. The spatial resolution (in meters) used to rasterize each slice.
#'
#' @return A list of three raster objects:
#' \describe{
#'   \item{\code{RasterMetric}}{A list of binary rasters where cells with at least one point are set to 1, others are NA.}
#'   \item{\code{RasterMetricSum}}{A raster layer representing the number of slices where each cell contains at least one point (sum of binary rasters).}
#'   \item{\code{RasterMetricCount}}{A raster layer representing the total number of points per cell summed across all slices.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Example usage with LAS data
#'   library(lidR)
#'   las <- readLAS("path_to_las_file.las")
#'   slices <- Slice3DByHeight(las, Z_mininum = 1, Z_maximum = 10, voxel_size = 0.2)
#'   raster_metrics <- Rasterize3DSlices(slices[[1]], voxel_size = 0.2)
#'   plot(raster_metrics[[2]])  # Raster showing the number of active layers per cell
#' }
#'
#' @import foreach
#' @import raster
#' @export
Rasterize3DSlices <- function(LasData,voxel_size) {
  Rasterize3DSlicesArch(LasData,voxel_size)
}

#' Identify Potential Tree Stem Locations from 3D Raster Metrics
#'
#' This function identifies potential tree stem locations based on multi-layer voxel raster data
#' derived from LiDAR point clouds. It evaluates stacked raster metrics representing height slices
#' and applies spatial and contextual filters to locate the most likely positions of tree stems.
#' This function is typically used after voxelizing point clouds with a consistent voxel size.
#'
#' @param xyz_data A data frame containing 3D point coordinates with columns `x`, `y`, and `z`.
#' @param RasterMetric A list of binary raster layers (one per voxel slice), indicating presence of points in each voxel layer.
#' @param RasterMetricSum A raster layer representing the number of non-empty slices across all voxel layers per cell.
#' @param RasterMetricCount A raster layer representing the sum of all values (e.g., point counts or intensities) across voxel layers per cell.
#' @param voxel_size A numeric value indicating the voxel resolution used to generate the raster layers.
#'
#' @return A data frame containing the filtered `xyz_data` rows that are classified as potential tree stem positions.
#'
#' @examples
#' \dontrun{
#' # Assuming voxel-based rasters have been calculated previously
#' potential_stems <- IdentifyPotentialTreeLocations(
#'   xyz_data = point_cloud_xyz,
#'   RasterMetric = voxel_layers,
#'   RasterMetricSum = raster_metric_sum,
#'   RasterMetricCount = raster_metric_count,
#'   voxel_size = 0.5
#' )
#' }
#'
#' @export
IdentifyPotentialTreeLocations <- function(xyz_data, RasterMetric, RasterMetricSum, RasterMetricCount, voxel_size) {
  IdentifyPotentialTreeLocationsArch(xyz_data, RasterMetric, RasterMetricSum, RasterMetricCount, voxel_size)
}

#' Tree Detection and Diameter Measurement from 3D Point Cloud Data
#'
#' This function detects tree trunks and estimates their diameters using a combination of DBSCAN clustering and RANSAC-based circle fitting. It is designed for analyzing 3D point cloud data from forest environments and extracting tree-level structural metrics such as DBH and basal area.
#' The function performs density-based clustering, filters clusters by minimum point count, fits circles using RANSAC, and returns a shapefile containing estimated diameters and basal areas for each detected tree.
#'
#' @param data_ref A data frame or matrix with 3D point cloud coordinates (x, y, z) representing tree points.
#' @param eps_value The neighborhood radius parameter (epsilon) for the DBSCAN clustering algorithm.
#' @param minPts_value The minimum number of points required to form a dense region in DBSCAN.
#' @param d1_m_minimum Minimum accepted diameter (in meters) for valid RANSAC-detected circles (used for filtering).
#' @param d1_m_maximum Maximum accepted diameter (in meters) for valid RANSAC-detected circles (used for filtering).
#' @param n_ransac_par The minimum number of points required to attempt RANSAC fitting on a cluster.
#' @param k_ransac_par Number of iterations RANSAC should use to select minimal point sets.
#' @param t_ransac_par Distance threshold for inlier classification in RANSAC fitting.
#' @param d_ransac_par Maximum distance from the model to be considered as an inlier (tolerance).
#'
#' @return A `sf` object (simple features) containing the detected trees with the following columns:
#' \item{id}{Unique ID of the detected cluster/tree.}
#' \item{dbh_cm}{Estimated diameter at breast height (in centimeters).}
#' \item{g_m2}{Estimated basal area (in square meters).}
#' \item{geometry}{Polygon geometries representing buffered circles around tree stem positions.}
#'
#' @examples
#' \dontrun{
#' library(dbscan)
#' library(dplyr)
#' library(sf)
#' # Simulate a small dataset
#' set.seed(1)
#' data_sim <- data.frame(x = runif(100, 0, 10),
#'                        y = runif(100, 0, 10),
#'                        z = runif(100, 0, 2))
#' # Run the detection and measurement function
#' tree_results <- DetectAndMeasureTrees(data_ref = data_sim,
#'                                       eps_value = 0.5,
#'                                       minPts_value = 10,
#'                                       d1_m_minimum = 0.05,
#'                                       d1_m_maximum = 1.0,
#'                                       n_ransac_par = 30,
#'                                       k_ransac_par = 100,
#'                                       t_ransac_par = 0.01,
#'                                       d_ransac_par = 0.02)
#' print(tree_results)
#' }
#'
#' @export
DetectAndMeasureTrees <- function(data_ref,eps_value,minPts_value,d1_m_minimum,
                                  d1_m_maximum, n_ransac_par, k_ransac_par,
                                  t_ransac_par, d_ransac_par) {
  DetectAndMeasureTreesArch(data_ref,eps_value,minPts_value,d1_m_minimum,
                            d1_m_maximum, n_ransac_par, k_ransac_par,
                            t_ransac_par, d_ransac_par)
}

#' Forest Stand Metrics from Tree Data
#'
#' This function calculates basic stand-level forestry metrics including total volume,
#' number of trees, and derived values per hectare. It takes as input a spatial object
#' containing detected trees with DBH measurements, and a polygon representing the
#' plot boundary.
#'
#' The tree volume is estimated using a generic allometric equation based on DBH:
#' \eqn{V = a \cdot \mathrm{DBH}^b} where \eqn{a = 0.0001078} and \eqn{b = 2.56}.
#'
#' @param result_final_base An `sf` object with at least a `dbh_cm` column representing the diameter at breast height (in centimeters) of individual trees.
#' @param tls_slice3 An `sf` polygon representing the plot area used to scale the metrics to per-hectare values.
#'
#' @return A list with the following components:
#' \item{total_volume_m3_plot}{Total tree volume (m³) in the input plot.}
#' \item{total_tree_n_plot}{Number of detected trees in the plot.}
#' \item{area_m2_plot}{Area of the plot in square meters.}
#' \item{total_volume_m3_ha}{Total tree volume (m³) scaled per hectare.}
#' \item{total_tree_n_ha}{Number of trees scaled per hectare.}
#' \item{area_m2_ha}{Area of the plot in hectares.}
#'
#' @examples
#' \dontrun{
#' # Example usage
#' library(sf)
#' poly_coords <- rbind(
#'   c(0,0), c(0,10), 
#'   c(10,10), c(10,0),
#'   c(0,0)
#' )
#' poly <- st_polygon(list(poly_coords))
#' tls_slice3 <- st_as_sf(
#'   st_sfc(poly), 
#'   crs = 3035
#' )
#'
#' @export
forestry_equations <- function(result_final_base, tls_slice3) {
  
  # 1. Area of the point cloud in hectares
  if (!inherits(result_final_base, "sf")) stop("point_cloud must be an 'sf' object")
  area_ha <- as.numeric(st_area(tls_slice3)) / 10000
  area_plot<- as.numeric(st_area(tls_slice3))
  # 2. Trees per hectare
  n_trees <- nrow(result_final_base)
  trees_per_ha <- n_trees / area_ha
  
  # 3. Generic volume function based on DBH
  volume_function <- function(dbh_cm) {
    a <- 0.0001078
    b <- 2.56
    return(a * (dbh_cm^b))
  }
  
  
  # 4. Calculate individual tree volume and total volume per hectare
  result_final_base$volume_tree_m3 <- volume_function(result_final_base[["dbh_cm"]])
  total_volume <- sum(result_final_base$volume_tree_m3, na.rm = TRUE)
  volume_per_ha <- total_volume / area_ha
  
  return(list(
    total_volume_m3_plot = round(total_volume,2),
    total_tree_n_plot = n_trees,
    area_m2_plot = round(area_plot,2),
    total_volume_m3_ha = round(volume_per_ha,2),
    total_tree_n_ha = round(trees_per_ha,2),
    area_m2_ha = round(area_ha,2)
  ))
}




