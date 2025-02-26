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

library(future)
start4 <- Sys.time()
# Set up parallel processing with ## workers
plan(multisession, workers = 8)

# Run the beta.fd.multidim() function in parallel
beta_FD <- future({
  beta.fd.multidim(
    sp_faxes_coord = sp_coords1[, paste("PC", 1:4, sep = "")],
    asb_sp_occ = as.matrix(presence_matrix),
    beta_family = c("Jaccard"),
    check_input = FALSE,
    details_returned = TRUE
  )
})

# Get the result
beta_FD_result <- value(beta_FD)

fin4 <- Sys.time()
Time4 <-  fin4-start4

saveRDS(Time4,"/scratch/user/llew0024/Time4.rds")
saveRDS(beta_FD_4PC, "/scratch/user/llew0024/beta_4PC")

# Set up parallel processing with 30 workers and 5 PCoAs
start5 <- Sys.time()
plan(multisession, workers = 30)

# Run the parallel computation
beta_FD_5PC <- future_map(1:30, function(i) {
  with_seed(123 + i, {
    beta.fd.multidim(
      sp_faxes_coord = sp_coords1[, paste("PC", 1:5, sep = "")],
      asb_sp_occ = as.matrix(presence_matrix),
      beta_family = c("Jaccard"),
      check_input = FALSE,
      details_returned = TRUE
    )
  })
})

fin5 <- Sys.time()
Time5 <- fin5-start5
saveRDS(Time5,"/scratch/user/llew0024/Time5.rds")
saveRDS(beta_FD_5PC, "/scratch/user/llew0024/beta_5PC")
  