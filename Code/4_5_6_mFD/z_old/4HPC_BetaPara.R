##Use mFD to apply PCoA and get some functional richness metrics, including all species diversity

# Add your custom library path
.libPaths(c("/home/llew0024/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()))

#Load packages
library(mFD)
library(dplyr)
library(tidyr)

#get the fish data
fish <- readRDS("/scratch/user/llew0024/fish.rds")
sp_coords1 <- readRDS("/scratch/user/llew0024/sp_coords1.rds")
presence_matrix <- readRDS("/scratch/user/llew0024/presence_matrix.rds")

# Set up parallel processing with 28 cores and 6 PCoAs
start6 <- Sys.time()

# Run the parallel computation
beta_FD_6PC <- beta.fd.multidim(
  sp_faxes_coord = sp_coords1[, paste("PC", 1:6, sep = "")],
  asb_sp_occ = as.matrix(presence_matrix),
  beta_family = c("Jaccard"),
  check_input = FALSE,
  details_returned = TRUE,
  betapart_para = TRUE,              # Enable parallelization
  betapart_para_opt = list(
    nc = 28,                          # Use 28 cores because there are 28 comparisons
    type = "PSOCK",                  # Use PSOCK clusters
    LB = TRUE,                       # Enable load-balancing
    size = 1                         # Number of tasks per worker
  )
)

fin6 <- Sys.time()
Time6 <- fin6-start6
saveRDS(Time6,"/scratch/user/llew0024/Time6.rds")
saveRDS(beta_FD_6PC, "/scratch/user/llew0024/beta_6PC")
  