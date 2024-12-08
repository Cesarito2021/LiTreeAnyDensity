# LiTreeAnyDensity
An R package for detecting and measuring tree attributes across 3D point clouds of any density

# ************************************************
#  Tutte le librerie
# ************************************************
suppressMessages(library(lidR))
suppressMessages(library(dplyr))
suppressMessages(library(raster))
suppressMessages(library(sp))
suppressMessages(library(terra))
suppressMessages(library(sf))
suppressMessages(library(dbscan))
suppressMessages(library(foreach))
suppressMessages(library(purrr))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library(microbenchmark))

# Set number of cores to use
numCores <- detectCores() - 1  # Usa tutti i core disponibili meno uno
cl <- makeCluster(numCores)
registerDoParallel(cl)

# ************************************************
#  Source Codes
# ************************************************
source("C:/Users/calvi/OneDrive/Desktop/Iphone_metodology/Rcode/baseScripts.R")

# ************************************************
#  library
# ************************************************
setwd("C:/Users/calvi/OneDrive/Desktop/Iphone_metodology/Dati")

# nome input
nome_sito = "iphone_M3_norm"
nome_sito_slice  = "iphone_M3_norm_slice"
# nome_sito_output = "parte_CC"

Z_minimum = 1
Z_maximum = 1.5
voxel_size = 0.01

# Parameters for DBSCAN and RANSAC
eps_value     = 0.01  # DBSCAN eps
minPts_value  = 3     # DBSCAN minPts
d1_m_minimum  = 0.07  # RANSAC min diameter
d1_m_maximum  = 0.9   # RANSAC max diameter
n_ransac_par = 20     # RANSAC iterations
k_ransac_par = 10     # RANSAC k parameter
t_ransac_par = 0.7    # RANSAC t parameter
d_ransac_par = 10     # RANSAC d parameter
crs_ransac_par = 32633  # CRS for RANSAC

# ************************************************
#  Import TLS data
# ************************************************
tls_slice  <- readTLSLAS(paste0(nome_sito_slice, ".las"))
tls  <- readTLSLAS(paste0(nome_sito, ".las"))

# ******************************************************
#           Slicing the point clouds in horizontal
# ******************************************************
start.time_step1 <- Sys.time()
result_step1 = Slice3DByHeight(tls_slice, Z_minimum, Z_maximum, voxel_size)
end.time_step1 <- Sys.time()
time.taken_step1 <- end.time_step1 - start.time_step1
time.taken_step1 

# ******************************************************
#           Matematical analysis in each Horizontal slice
# ******************************************************
start.time_step2 <- Sys.time()
result_step2 = Rasterize3DSlices(result_step1[[1]], voxel_size)
end.time_step2 <- Sys.time()
time.taken_step2 <- end.time_step2 - start.time_step2
time.taken_step2

# ******************************************************
#           Selecting points falling in trunks cylinders
# ******************************************************
start.time_step3 <- Sys.time()
result_step3_improved = IdentifyPotentialTreeLocations(result_step1[[2]], result_step2[[1]], result_step2[[2]], result_step2[[3]], voxel_size)
end.time_step3 <- Sys.time()
time.taken_step3 <- end.time_step3 - start.time_step3
time.taken_step3 

# **************************************************************************
#           Detecting and measuring circles of trunks and tree height
# **************************************************************************
start.time_step4 <- Sys.time()
result_final_base <- DetectAndMeasureTrees(result_step3_improved, eps_value, minPts_value, d1_m_minimum, d1_m_maximum,
                                           paste(nome_sito_output), crs_ransac_par, n_ransac_par, k_ransac_par,
                                           t_ransac_par, d_ransac_par, 530)
end.time_step4 <- Sys.time()
time.taken_step4 <- end.time_step4 - start.time_step4
time.taken_step4

# Plot for base data with customized circle colors and sizes
plot(result_final_base$centroid_x, result_final_base$centroid_y,
     main = "Base Data",
     xlab = "X Coordinate", ylab = "Y Coordinate",
     col = "red", pch = 16, cex = 0.8)  # red color, filled circles, smaller size

stopCluster(cl)
getwd()

# # Example of printing times with units
cat("SLICING time secs:")
time.taken_step1
cat("SLICING PROCESS time min:")
time.taken_step2
cat("SLICING CLEANING time min:")
time.taken_step3
cat("DETECTING/MEASURING TREES time min:")
time.taken_step4

# **************************************************************************
#                                THE END 
# **************************************************************************

