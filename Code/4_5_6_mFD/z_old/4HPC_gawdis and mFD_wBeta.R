##Use mFD to apply PCoA and get some functional richness metrics, including all species diversity

#Load packages
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(cowplot)

#get the fish data
fish <- readRDS("/scratch/user/llew0024/fish.rds")

#separate trait data, dropping SL because of inconsistency in how it is measured and mandible because it is highly skewed
trait_names <- c("Species","BodyShapeI","BodyShapeII","TL","HL","ED","POL","BD","PosofMouth","eye.position","spiracle","caudal.fin.shape") #"mandible"
trait <- fish[,trait_names]
trait <- unique(trait)
rownames(trait) <- trait$Species
trait$Species <- NULL
details <- fish[,names(fish)%in%c("Species","community", "habitat", "degsFromEquator")]

#fix BodyShapeI - fusiform had been split into two groups
trait$BodyShapeI <- as.character(trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("fusi",trait$BodyShapeI), "fusiform",trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("short",trait$BodyShapeI), "shortAndOrDeep",trait$BodyShapeI)

#combine eye.position front-facing and raised/top of head (only Bothriolepis canadensis has front-facing, but they are quite raised)
trait$eye.position <- as.character(trait$eye.position)
trait$eye.position <- ifelse(trait$eye.position=="front-facing","raised/top of head",trait$eye.position)

################################################################################
#check if the quantitative traits should be log-transformed
nums <- trait[sapply(trait, is.numeric)] 
#plot histograms
lapply(names(nums), function(col_name) {
  hist(nums[[col_name]], main = col_name)
})
#make the transformations where needed
#TL
nums$TL <-  log(nums$TL)
names(nums)[names(nums)=="TL"] <- "logTL"
trait$TL <-  log(trait$TL)
names(trait)[names(trait)=="TL"] <- "logTL"
#POL
nums$POL <-  log(nums$POL)
names(nums)[names(nums)=="POL"] <- "logPOL"
trait$POL <-  log(trait$POL)
names(trait)[names(trait)=="POL"] <- "logPOL"
#check distribution again
lapply(names(nums), function(col_name) {
  hist(nums[[col_name]], main = col_name)
})

#make assemblage matrix
# Add a column for presence (1) for each species in each community
details <- details %>% mutate(presence = 1)

# Transform the data frame into a wide format
presence_matrix <- details %>%
  spread(key = Species, value = presence, fill = 0)

# remove the community column and community details from the matrix
row.names(presence_matrix) <- presence_matrix$community
cd <- presence_matrix[,1:2]
presence_matrix <- presence_matrix[, !names(presence_matrix)%in%c("community","habitat","degsFromEquator")]

#need trait detail data frame first, with character columns changed to factor and row names = species
trait_det <- data.frame(trait_name=names(trait) ,trait_type=ifelse(sapply(trait,is.numeric),"Q","N"),trait_weight = 1, fuzzy_name = NA)
traits <- lapply(trait, function(x) if(class(x) == "character") as.factor(x) else x)
traits <- as.data.frame(traits)
rownames(traits) <- rownames(trait)
#check content of the data frames
sp.tr.summary(sp_tr = traits, tr_cat = trait_det)

##get distances using gawdis###################################################################
GD_all <- gawdis(traits) #no issue with negative weights, but have unbalanced distributions involving some traits
attr(GD_all,"correls") #contribution are equal
attr(GD_all,"weights") #all have positive weights

#  #figure out which traits to keep - must have no negative weights (but unbalanced distribution is okay - indicates distribution of values that is heavily skewed or that these traits have many identical values for most of the species - it means that rare categories only have weak contriutions to distances#
#  num_traits <- ncol(traits)
#  max_combinations <- vector("list", num_traits)
#  max_no_warning <- 0
#  max_combination_names <- list()
# 
# # Function to apply gawdis and catch warnings
# safe_gawdis <- function(combination) {
#   warnings <- NULL
#   result <- withCallingHandlers(
#     gawdis(traits[, combination, drop = FALSE]),
#     warning = function(w) warnings <<- c(warnings, w$message)
#   )
#   list(result = result, warnings = warnings)
# }
# 
# # Initialize max_combinations for each possible length
# for (i in 8:ncol(traits)) {
#   max_combinations[[i]] <- list()
# }
# 
# # Loop over all combinations starting from 8 - it takes a while
# for (i in 8:ncol(traits)) {
#   combinations <- combn(num_traits, i, simplify = FALSE)
# 
#   for (combination in combinations) {
#     result <- safe_gawdis(combination)
# 
#     if (is.null(result$warnings)) {  # No warnings
#       combination_length <- length(combination)
#       if (combination_length > max_no_warning) {
#         max_no_warning <- combination_length
#         max_combination_names <- list(colnames(trait)[combination])
#       } else if (combination_length == max_no_warning) {
#         max_combination_names <- c(max_combination_names, list(colnames(trait)[combination]))
#       }
#     }
#   }
# }
# 
# # max_combination_names contains the names of traits in the models with the most traits that didn't produce warnings
# max_combination_names

#compare combinations that don't have unbalanced data 
#the combinations:
c1 <- c("BodyShapeI","BodyShapeII","logTL","HL","ED","logPOL","BD","PosofMouth","caudal.fin.shape")
c2 <- c("BodyShapeI","logTL","HL","ED","logPOL","BD","PosofMouth","eye.position", "caudal.fin.shape")
c3 <- c("BodyShapeI", "logTL","HL","ED","logPOL","BD","PosofMouth","spiracle","caudal.fin.shape")
combos <- list(c1,c2,c3)

#############################################################################
#make list of trait data sets with the different trait combination
traits_cut <- list()
for(i in 1:length(combos)){
  traits_cut[[i]] <- traits[,combos[[i]]]}
  
#now get the distances 
gd <- lapply(traits_cut,gawdis) #using gawdis function directly
CandW <- function(x) {
  conts <-attr(x,"correls") #contribution are equal
  wghts <- attr(x,"weights") # weight automatically adjusted so traits have same contribution
  results <- list(conts,wghts)
  return(results)}
check <- lapply(gd,CandW)
check_all <- CandW(GD_all)

#Compute multimensional functional spaces (PCoA) and assess their quality
qual <- lapply(gd, quality.fspaces,fdendro = "average",maxdim_pcoa = 487,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE))
qual1 <- quality.fspaces(sp_dist = GD_all, fdendro = "average",maxdim_pcoa = 487,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 

#check if any species have 0 distance
duplicates_all <- lapply( traits_cut, function(x){duplicated(x) | duplicated(x, fromLast = TRUE)})
table(duplicates_all)
lapply(duplicates_all, table) #no duplicates in terms of trait sets

#check mad index to identify best functional space
# retrieve the functional space associated with minimal quality metric: 
lapply(qual,function(x){apply(x$quality_fspaces, 2, which.min)}) #different measures suggest different number of dimensions (6 or 7)
apply(qual1$quality_fspaces,2,which.min)

#plot it to see
library("magrittr")

lapply(qual, function(x){
x$"quality_fspaces" %>%
  tibble::as_tibble(rownames = "Funct.space") %>%
  tidyr::pivot_longer(cols =! Funct.space, names_to = "quality_metric", values_to = "Quality") %>%
  ggplot2::ggplot(ggplot2::aes(x = Funct.space, y = Quality, 
                               color = quality_metric, shape = quality_metric)) +
  ggplot2::geom_point()}) 
#5 to 9 PCoAs looks good

#Alternatively, choose by how much variation explained by each PCoA
lapply(qual, function(x){eigs <- x$details_fspaces$pc_eigenvalues
Esum <- sum(eigs$Eigenvalues)
eigs$perc <- eigs$Eigenvalues/Esum
cbind(rownames(eigs),cumsum(eigs$perc))
plot(rownames(eigs),cumsum(eigs$perc))
print(cumsum(eigs$perc)[7])}) #7 eigs  explain > 72%
eigs <- qual1$details_fspaces$pc_eigenvalues
Esum <- sum(eigs$Eigenvalues)
eigs$perc <- eigs$Eigenvalues/Esum
cbind(rownames(eigs),cumsum(eigs$perc))
plot(rownames(eigs),cumsum(eigs$perc))
print(cumsum(eigs$perc)[8]) #8 eigs explain 73% when all traits are included

#position of species
sp_coords <- lapply(qual, function(x){x$details_fspaces$sp_pc_coord})
sp_coords1 <- qual1$details_fspaces$sp_pc_coord

#make file for exporting and using to calculate overlap
species_coords <- lapply(sp_coords, function(x) {
  spc <- data.frame(x)
  pm <- data.frame(t(presence_matrix))
  spc$species <- row.names(spc)
  pm$species <- row.names(pm)
  pm$species <-  gsub("\\.", " ", pm$species)
  spc <- merge(spc,pm, by=("species"), all.x=TRUE)})
#saveRDS(species_coords, "Data/PCoAs/species_coordinates.RDS")

#see correlations between traits and axes #can handle up to 10 traits, produces data frames and plots
Tcorrs <- mapply(function(x,y){
  obs<-traits.faxes.cor(sp_tr = x, sp_faxes_coord = y[,paste("PC",1:7, sep="")], plot = TRUE)}, 
  x=traits_cut, y=sp_coords)
           
#get correlations between traits and axes 
Tcorrs1 <- traits.faxes.cor(sp_tr = traits, sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")])
Tcorrs1 <- Tcorrs1[order(Tcorrs1$axis, Tcorrs1$value), ]
#save it
#write.csv(Tcorrs1,"~/Dropbox/Devonian fish/manuscript/ProcB revision/tables/PC_trait_correls.csv", row.names = FALSE)

#calculate functional diversity indices
#need matrix of 1s and 0s where row = community and column = fish species -> the presence_matrix
#alpha_FD <- lapply(sp_coords, function(x){alpha.fd.multidim(sp_faxes_coord = x[,paste("PC",1:7, sep="")], asb_sp_w = as.matrix(presence_matrix))}) #can scale (scaling = TRUE) so values are between 0 and 1, but different indices are squashed into different portions of this range; can instead set mean to 0 and sd to 1
#alpha_FD3 <- alpha.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")], asb_sp_w = as.matrix(presence_matrix))

#calculate Beta diversity m- personal computer can't handle it with that many demensions. Run on HPC
#beta_FD <- lapply(sp_coords, function(x){beta.fd.multidim(sp_faxes_coord = x[,paste("PC",1:7, sep="")], asb_sp_occ = as.matrix(presence_matrix), beta_famil = c("Jaccard"), check_input = TRUE, details_returned = TRUE)}) #can scale (scaling = TRUE) so values are between 0 and 1, but different indices are squashed into different portions of this range; can instead set mean to 0 and sd to 1
#beta_FD3 <- beta.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:4, sep="")], asb_sp_occ = as.matrix(presence_matrix),beta_famil = c("Jaccard"), check_input = TRUE, details_returned = TRUE)

library(furrr)
library(withr)
start4 <- Sys.time()
# Set up parallel processing with 30 workers
plan(multisession, workers = 30)

# Run the parallel computation
beta_FD_4PC <- future_map(1:30, function(i) {
  with_seed(123 + i, {
    beta.fd.multidim(
      sp_faxes_coord = sp_coords1[, paste("PC", 1:4, sep = "")],
      asb_sp_occ = as.matrix(presence_matrix),
      beta_family = c("Jaccard"),
      check_input = FALSE,
      details_returned = TRUE
    )
  })
})

fin4 <- Sys.time()
Time4 <- start4-fin4
saveRDS(Time4,"Time4.rds")

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
Time5 <- start5-fin5
saveRDS(Time5,"Time5.rds")

saveRDS(beta_FD_5PC, "/scratch/user/llew0024/beta_5PC")
  