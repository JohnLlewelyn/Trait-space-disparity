#betapart (Baselga et al.)

library(betapart)

setwd("/Users/llew0024/Dropbox/Devonian fish/manuscript/ProcB revision/GitHub_code&data_revised")
fish <- readRDS("Data/HPC/fish.rds")
sp_coords1 <- readRDS("Data/HPC/sp_coords1.rds")
presence_matrix <- readRDS("Data/HPC/presence_matrix.rds")

#fix object
sp_coords1 <- data.frame(sp_coords1)
sp_coords1 <- sp_coords1[names(presence_matrix),]
sp_coords1 <- sp_coords1[,1:8]

#run it
func_betapart_obj <- functional.betapart.core(
  x = presence_matrix,
  traits = sp_coords1,
  return.details = TRUE
)
