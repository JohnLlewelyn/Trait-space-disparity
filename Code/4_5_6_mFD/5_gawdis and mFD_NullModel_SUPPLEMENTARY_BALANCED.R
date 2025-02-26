#Use mFD to apply PCoA and get some functional richness metrics, controlling for species using standardized effect sizes
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(grid)
library(gridExtra)
library(cowplot)

#Adjust file path at"###"

#get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source("###/Code/2_3_tidy data/2_Devonian_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

setwd(wkD)

#get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source("###/Code/2_3_tidy data/3_Devonian_tidy_Gogo_Miguasha.R", local = new_env)
dv <- get("mg", envir = new_env)
rm(new_env)

#remove the extra stuff
rm(list= ls()[! (ls() %in% c('fish','dv'))])

#stick it together
fish <- fish[,names(dv)]
fish <- rbind(fish,dv)

#remove Little Rock Lake; too few fish
fish <- fish[fish$community!="traits_Little_Rock_Lake",]

#separate trait data, dropping SL because of inconsistency in how it is measured and mandible because it is highly skewed
trait_names <- c("Species","BodyShapeI","BodyShapeII","TL","HL","ED","POL","BD","PosofMouth","eye.position","spiracle","caudal.fin.shape" )
trait <- fish[,trait_names]
trait <- unique(trait)
rownames(trait) <- trait$Species
trait$Species <- NULL
details <- fish[,names(fish)%in%c("Species","community", "habitat", "degsFromEquator")]

#fix BodyShapeI - fusiform had been split into two groups
trait$BodyShapeI <- as.character(trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("fusi",trait$BodyShapeI), "fusiform",trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("short",trait$BodyShapeI), "shortAndOrDeep",trait$BodyShapeI)

#combine eye.position front-facing and raised/top of head (only Bothriolepis canadensis has front-facing, and they are quite raised)
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
GD_all <- gawdis(traits) #no issue with negative weights, but still have unbalanced distributions for some traits
attr(GD_all,"correls") #contribution are equal
attr(GD_all,"weights") 

#Compute multimensional functional spaces (PCoA) and assess their quality
qual1 <- quality.fspaces(sp_dist = GD_all, fdendro = "average",maxdim_pcoa = 487,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 

#compare combinations that don't have unbalanced data #see 4_gawdis and mFD_allSP.R for how these were identified
#the combinations:
c1 <- c("BodyShapeI","BodyShapeII","logTL","HL","ED","logPOL","BD","PosofMouth","caudal.fin.shape")
c2 <- c("BodyShapeI","logTL","HL","ED","logPOL","BD","PosofMouth","eye.position", "caudal.fin.shape")
c3 <- c("BodyShapeI", "logTL","HL","ED","logPOL","BD","PosofMouth","spiracle","caudal.fin.shape")
combos <- list(c1,c2,c3)

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
qual <- lapply(gd, quality.fspaces,fdendro = "average",maxdim_pcoa = 487,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE))

#position of species
sp_coords1 <- qual1$details_fspaces$sp_pc_coord
sp_coords <- lapply(qual, function(x){x$details_fspaces$sp_pc_coord})

################################################################################
#observed values
obs_alpha <- alpha.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")], asb_sp_w = as.matrix(presence_matrix))
obs_alphaL <- lapply(sp_coords, function(x) alpha.fd.multidim(sp_faxes_coord = x[,paste("PC",1:8, sep="")], asb_sp_w = as.matrix(presence_matrix)))

#Generate null model distributions
#Set parameters for null model
n_simulations <- 100  # Number of null simulations

#function for shuffling matrix, keeping row and column sum unchanged
# Function to shuffle a binary matrix while strictly preserving row and column sums
shuffle_matrix_fixed <- function(mat, n_swaps) {
  mat <- as.matrix(mat)  # Ensure input is a matrix
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  for (swap in 1:n_swaps) {
    # Randomly select two rows
    rows <- sample(1:n_rows, 2, replace = FALSE)
    # Randomly select two columns
    cols <- sample(1:n_cols, 2, replace = FALSE)
    
    # Extract the 2x2 submatrix
    submatrix <- mat[rows, cols]
    
    # Check if the submatrix is swappable (it must contain exactly two 1s and two 0s)
    if (sum(submatrix) == 2 && all(rowSums(submatrix) == 1) && all(colSums(submatrix) == 1)) {
      # Swap the 1s and 0s
      mat[rows[1], cols[1]] <- 1 - mat[rows[1], cols[1]]
      mat[rows[1], cols[2]] <- 1 - mat[rows[1], cols[2]]
      mat[rows[2], cols[1]] <- 1 - mat[rows[2], cols[1]]
      mat[rows[2], cols[2]] <- 1 - mat[rows[2], cols[2]]
    }
  }
  
  return(mat)
}
# ####SHUFFLING # Generate null model distributions separately for each set of species coordinates
###It takes a while, so results have also been saved and can be directly read in (below)
# # Set parameters for null model
# n_simulations <- 1000  # Number of null simulations
# n_swaps <- 10000      # Number of swaps during shuffling
# 
# # Initialize a list to store null results for each set of species coordinates
# null_fd_list <- list()
# 
# # Perform null model simulations for each set of species coordinates
# start_time <- Sys.time()
# set.seed(89)  # For reproducibility
# 
# for (coord_set in seq_along(sp_coords)) {
#   # Initialize matrix to store null FRic values for the current coordinate set
#   null_fd <- list()
#   
#   for (i in 1:n_simulations) {
#     # Shuffle species labels across the entire species pool
#     shuffled_asb_sp_w <- shuffle_matrix_fixed(as.matrix(presence_matrix), n_swaps = n_swaps)
#     
#     # Check row and column sums are maintained (optional debugging)
#     # print(rowSums(shuffled_asb_sp_w))
#     # print(table(colSums(shuffled_asb_sp_w)))
#     # print(shuffled_asb_sp_w[, 1:4])
#     
#     # Compute FRic for shuffled communities using the current set of species coordinates
#     fd <- mFD::alpha.fd.multidim(
#       sp_faxes_coord = sp_coords[[coord_set]][, paste("PC", 1:8, sep = "")],
#       asb_sp_w = shuffled_asb_sp_w
#     )
#     
#     # Store null FRic values
#     null_fd[[i]] <- fd[[1]][, 1:which(colnames(fd[[1]]) == "fspe")]
#     print(paste("Simulation", i, "for coordinate set", coord_set))
#   }
#   
#   # Save null results for the current coordinate set
#   null_fd_list[[coord_set]] <- null_fd
# }
# 
# # Save the null results for all coordinate sets
# saveRDS(null_fd_list, "Data/null_results_sp_coords.RDS")
# end_time <- Sys.time()
# end_time - start_time  # ~ 9 hours

null_fd <- readRDS("###/Data/data_for_plots/null_results_sp_coords.RDS")

#Calculate means and SDs across lists e.g., mead and SD for each communities Fric
#mean
mean_fd <- lapply(null_fd, function(x) Reduce("+", x) / length(x))

#SD
n_rows <- lapply(null_fd, function(x) nrow(x[[1]]))
n_cols <- lapply(null_fd, function (x) ncol(x[[1]]))

# Initialize an empty matrix for the standard deviations
sd_fd <- matrix(NA, nrow = n_rows[[1]], ncol = n_cols[[1]])
sd_fd <- list(sd_fd, sd_fd, sd_fd)

# Loop through rows and columns
for (i in 1:3){
for (j in seq_len(n_rows[[i]])) {
  for (k in seq_len(n_cols[[i]])) {
    # Extract the values for the same cell across all data frames
    cell_values <- sapply(null_fd[[i]], function(df) df[j, k])
    # Calculate standard deviation
    sd_fd[[i]][j, k] <- sd(cell_values)
  }
}
}

sd_fd <- lapply(sd_fd, data.frame)
names(sd_fd) <- names(obs_alpha[[1]])[1:9]

# Get the new column names from obs_alpha[[1]]
new_column_names <- names(obs_alpha[[1]])[1:9]

# Rename the columns in each data frame in the list
sd_fd <- lapply(sd_fd, function(df) {
  colnames(df)[1:9] <- new_column_names  # Assign new column names to the first 9 columns
  return(df)
})

wkD <- "###/"
saveRDS(obs_alphaL, paste(wkD,"/data/data_for_plots/Plot_orig_combos.RDS", sep = ""))

# Calculate effect size
ES <- matrix(NA, nrow = 8, ncol = 8)
ES <- list(ES,ES,ES)
for (i in 1:3){
ES[[i]] <- obs_alphaL[[i]][[1]][,which(names(obs_alphaL[[i]][[1]])=="fdis"):which(names(obs_alphaL[[i]][[1]])=="fspe")] - mean_fd[[i]][,which(names(mean_fd[[i]])=="fdis"):which(names(mean_fd[[i]])=="fspe")]
}

# Divide by standard deviation to standardise
SES <- matrix(NA, nrow = 8, ncol = 8)
SES <- list(SES,SES,SES)
for (i in 1:3){
SES[[i]] <- ES[[i]]/sd_fd[[i]][,which(names(sd_fd[[i]])=="fdis"):which(names(sd_fd[[i]])=="fspe")]
}  

#add site
SES <- lapply(SES, function(x) data.frame(x, site =rownames(x)))

# Convert the data frame from wide to long format
longFM <- function(x) {
  x$habitat <- ifelse(x$site%in%c("Gogo","traits_Caribbean","traits_Chile_reef"), "reef",
                      ifelse(x$site%in%c("Miguasha","traits_Ythan","Santa_Cruz_Channel"), "estuary","fresh water"))
  x$group <- ifelse(x$site%in%c("Miguasha", "Gogo"), "Devonian",
                    ifelse(x$site%in%c("Santa_Cruz_Channel","traits_BracoMorto","traits_Caribbean"),"tropical","temperate/sub-trop"))
  long_indsC <- pivot_longer(x,
                             cols = c("fdis","feve","fric","fmpd","fnnd","fdiv","fori","fspe"),
                             names_to = "metric",
                             values_to = "value")
}

# long_indsct <- lapply(indsct, longFM)
long_indsC1 <- lapply(SES,longFM)

#tidy site names
Fnames <- function(x){
  x$site <- gsub("traits_","",x$site)
  x <- as.data.frame(x)
  add1 <- x[x$group%in%"Devonian",] #to make Devonian sites stand out
  x <- rbind(x,add1)
  return(x)}

#long_indsct <- lapply(long_indsct,Fnames)
long_indsC1 <- lapply(long_indsC1,Fnames)

#Fix names of metrics
#functional dispersion, functional evenness, functional richness
Mnames <- function(x) {
  x$metric <- ifelse(x$metric=="fdis","dispersion",
                     ifelse(x$metric=="feve", "eveness",
                            ifelse(x$metric=="fric", "richness",
                                   ifelse(x$metric=="fmpd", "pairwise\ndist.",
                                          ifelse(x$metric=="fnnd", "nearest\nneighbour",
                                                 ifelse(x$metric=="fdiv", "divergence",
                                                        ifelse(x$metric=="fori", "originality","specialization")))))))
  return(data.frame(x))}

# 
long_indsC1 <- lapply(long_indsC1, Mnames)

#Create plots in ggplot
#colour indicates ancient , modern, tropical, temperate
#shape indicates habitat type
#function for plotting the different metrics
plotF <- function(long_indsC1) {
  ggplot(long_indsC1, aes(x = metric, y = value, shape = habitat, color = group)) +
    geom_point(size = 5, alpha = 0.6) +  # Your existing points
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(angle = 0, size = 18),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      plot.title = element_text(size = 18)) +
    guides(size = FALSE, alpha = FALSE) +
    labs(x = "functional diversity metric",
         y = "metric value")
}

# combo_plots <- lapply(long_indsct, plotF)
metricsPlots1 <- lapply(long_indsC1, plotF)

#all traits figure
pdf("###/Allmetrics_Combos_SES.pdf",width=34, height=10)
 final_plot <- plot_grid(
   metricsPlots1[[1]] + theme(legend.position = "none"),
   metricsPlots1[[2]] + theme(legend.position = "none"),
   metricsPlots1[[3]],
   ncol = 3, nrow = 1, rel_widths = c(3, 3, 3.7)  # Specify the relative widths
 )
 print(final_plot)
dev.off()

# #cut to five focal metrics
metrics <- c("richness","nearest\nneighbour","eveness","divergence","specialization")
lc1 <- lapply(long_indsC1, function(x) x <- x[x$metric%in%metrics,])

lc1 <- lapply(lc1, function(df) {
  df$metric <- factor(df$metric, levels = c("richness", "nearest\nneighbour", "eveness", "divergence", "specialization"))
  return(df)
})

Plot1 <- lapply(lc1, plotF)

#combine the combos in a single plot
pdf("###/Sfigure 4_5metrics_Combos_controlDiversity.pdf",width=28, height=10)
final_plot <- plot_grid(
  Plot1[[1]] + theme(legend.position = "none"),
  Plot1[[2]] + theme(legend.position = "none"),
  Plot1[[3]],
  ncol = 3, nrow = 1, rel_widths = c(3, 3, 3.9)  # Specify the relative widths
)
print(final_plot)
dev.off()

saveRDS(Plot1,"###/Data/data_for_plots/Plot1_SESsupp.RDS")
