#' @importFrom purrr map map_dbl map_df
#' @importFrom Rvcg vcgClean
#' @importFrom dbscan dbscan
#' @importFrom lidR readLAS filter_poi
#' @importFrom openxlsx write.xlsx
#' @importFrom parallel detectCores
#' @importFrom terra rast
#' @importFrom dplyr %>% filter group_by mutate rename bind_rows left_join n ungroup group_split
#' @importFrom sf st_area st_as_sf st_sfc st_polygon
#' @importFrom stats optim prcomp
#' @importFrom lidR LAS filter_poi
NULL

utils::globalVariables(c(
  "x", "y", "z", "ID_Cluster", "id", "j", "k", 
  "d1_m", "dbh_cm", "g_m2", "geometry", "stem_cl",
  "i", "Z", "interval", "X", "Y"
))
if (!requireNamespace("purrr", quietly = TRUE)) {
  stop("Package 'purrr' needed but not installed.")
}

purrr::map(...)
