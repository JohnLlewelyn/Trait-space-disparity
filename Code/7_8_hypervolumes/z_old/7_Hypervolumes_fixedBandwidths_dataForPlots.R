#build hypervolumes and get metrics

library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(cowplot)
library(hypervolume)
library(colorBlindness)

#Adjust file path at"###/"
wkD <- "~/Dropbox/Devonian fish/manuscript/ProcB revision/GitHub_code&data_revised"
setwd(wkD)

#get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source("Code/2_3_tidy data/2_Devonian_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

wkD <- "~/Dropbox/Devonian fish/manuscript/ProcB revision/GitHub_code&data_revised"
setwd(wkD)

#get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source("Code/2_3_tidy data/3_Devonian_tidy_Gogo_Miguasha.R", local = new_env)
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
attr(GD_all,"correls") #contributions are equal
attr(GD_all,"weights") #all have positive weights

#Compute mPCoA and assess quality#######################################
#Compute multimensional functional spaces (PCoA) 
qual <- quality.fspaces(sp_dist = GD_all, fdendro = "average",maxdim_pcoa = 487,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 
#position of species
sp_coords <- qual$details_fspaces$sp_pc_coord

#pull the rows out for each assemblage using presence_matrix and sp_coords or spco
#Initialize an empty list to store the data frames for each site
site_data_frames1 <- list()
# Loop through each site in the presence matrix
for (site in rownames(presence_matrix)) {
  # Identify species present at this site
  species_present <- colnames(presence_matrix)[presence_matrix[site, ] == 1]
  # Extract the rows for these species from spco1
  site_df <- sp_coords[species_present, ]
  # Add the new data frame to the list, named by the site
  site_data_frames1[[site]] <- site_df
}

#make a separate site_data for Devonian versus modern and a global set
site_data_DevVSmod <- list(Devonian = rbind(site_data_frames1$Gogo,site_data_frames1$Miguasha), 
                           Modern = rbind(site_data_frames1$Santa_Cruz_Channel,
                                 site_data_frames1$traits_BracoMorto, 
                                 site_data_frames1$traits_Caribbean,
                                 site_data_frames1$traits_Chile_reef,
                                 site_data_frames1$traits_Nepean,
                                 site_data_frames1$traits_Ythan))

site_data_global <- list(Global = do.call(rbind, site_data_frames1))
 
#get rid of duplicates
site_data_DevVSmod <- lapply(site_data_DevVSmod, unique)
site_data_global <- lapply(site_data_global, unique)

#Use first 8 PCoAs to build hypervolumes, get the average bandwidths for each axis and then set them as fixed bandwidths (so the same are used for each assemblage)
#PCoA axes are roughly normally distributed, so use Guassian method
set.seed(123) 
start_time <- Sys.time()
hvsprep <- lapply(site_data_frames1, function(x){x <-x[,1:8]})
hvs1 <- lapply(hvsprep, hypervolume_gaussian)
end_time <- Sys.time()
end_time - start_time #1.145 minutes

set.seed(321) 
start_time <- Sys.time()
DMprep <- lapply(site_data_DevVSmod, function(x){x <-x[,1:8]})
hvs2 <- lapply(DMprep, hypervolume_gaussian)
end_time <- Sys.time()
end_time - start_time

set.seed(567) 
start_time <- Sys.time()
Gprep <- lapply(site_data_global, function(x){x <-x[,1:8]})
hvsG <- lapply(Gprep, hypervolume_gaussian)
end_time <- Sys.time()
end_time - start_time

#get the bandwidths used for each community for each PCoA
bandw <- data.frame(site=character(), bw1=numeric(),bw2=numeric(),bw3=numeric(),bw4=numeric(),
                    bw5=numeric(),bw6=numeric(),bw7=numeric(),bw8=numeric(), stringsAsFactors=FALSE)
# Loop through each hypervolume in hvs
for(i in 1:length(hvs1)) {
  # Append a new row to rdf
  bandw <- rbind(bandw, data.frame(site=names(hvs1)[i], bw1=hvs1[[i]]@Parameters$kde.bandwidth[1],
                    bw2=hvs1[[i]]@Parameters$kde.bandwidth[2],bw3=hvs1[[i]]@Parameters$kde.bandwidth[3],
                    bw4=hvs1[[i]]@Parameters$kde.bandwidth[4],bw5=hvs1[[i]]@Parameters$kde.bandwidth[5],
                    bw6=hvs1[[i]]@Parameters$kde.bandwidth[6],bw7=hvs1[[i]]@Parameters$kde.bandwidth[7],
                    bw8=hvs1[[i]]@Parameters$kde.bandwidth[8],stringsAsFactors=FALSE))
}
mns <-colSums(bandw[,2:9])/8

#recalculate hvs for each community with specified bandwidths
set.seed(123)
start_time <- Sys.time()
hvsprep <- lapply(site_data_frames1, function(x){x <-x[,1:8]})
bandwidths_fixed <- lapply(hvsprep, function(x) {x <-estimate_bandwidth(x,method="fixed", value=mns)})
hvs <- list()
for(i in 1: length(hvsprep)){
  hvs[[i]] <- hypervolume_gaussian(hvsprep[[i]], kde.bandwidth=bandwidths_fixed[[i]])
}
end_time <- Sys.time()
end_time - start_time #takes a minute

names(hvs) <- names(hvsprep)
saveRDS(hvs, "Data/data_for_plots/hvsCommunity.rds")
saveRDS(hvsprep, "Data/data_for_plots/hvsprepCommunity.rds")

############
#recalculate hvs for modern versus Devonian comparison with specified bandwidths
set.seed(123)
start_time <- Sys.time()
DMprep <- lapply(site_data_DevVSmod, function(x){x <-x[,1:8]})
bandwidths_fixedDM <- lapply(DMprep, function(x) {x <-estimate_bandwidth(x,method="fixed", value=mns)})
hvsDM <- list()
for(i in 1: length(DMprep)){
  hvsDM[[i]] <- hypervolume_gaussian(DMprep[[i]], kde.bandwidth=bandwidths_fixedDM[[i]])
}
end_time <- Sys.time()
end_time - start_time #takes 40 sec

#assign names
names(hvsDM) <- names(DMprep)

saveRDS(hvsDM, "Data/data_for_plots/hvsDM.rds")
saveRDS(DMprep, "Data/data_for_plots/DMprep.rds")

############
#Global hypervolume
set.seed(529)
start_time <- Sys.time()
Gprep <- lapply(site_data_global, function(x){x <-x[,1:8]})
bandwidths_fixedG <- lapply(Gprep, function(x) {x <-estimate_bandwidth(x,method="fixed", value=mns)})
hvsG <- list()
for(i in 1: length(Gprep)){
  hvsG[[i]] <- hypervolume_gaussian(Gprep[[i]], kde.bandwidth=bandwidths_fixedG[[i]])
}
end_time <- Sys.time()
end_time - start_time #takes 35 sec


saveRDS(hvsG, "Data/data_for_plots/hvsGlobal.rds")
saveRDS(Gprep, "Data/data_for_plots/Gprep.rds")

############
#cut to 400,000 random points for each
hvs <- lapply(hvs, function(hv) {
  hv@RandomPoints <- hv@RandomPoints[sample(nrow(hv@RandomPoints), 400000), ]
  return(hv) 
})

hvsG[[1]]@RandomPoints <- hvsG[[1]]@RandomPoints[sample(nrow(hvsG[[1]]@RandomPoints), 400000), ]

#make lists to plot
Gogo <- list(Global = hvsG[[1]], 'Gogo Reef' = hvs[["Gogo"]])
Miguasha <- list(Global = hvsG[[1]], 'Miguasha Estuary' = hvs[["Miguasha"]])
SantaCruz <- list(Global = hvsG[[1]], 'Santa Cruz Estuary' = hvs[["Santa_Cruz_Channel"]])
BracoMorto <- list(Global = hvsG[[1]], 'Braço Morto Acima and Abaixo' = hvs[["traits_BracoMorto"]])
Caribbean <- list(Global = hvsG[[1]], 'Caribbean Reef' = hvs[["traits_Caribbean"]])
Chile <- list(Global = hvsG[[1]], 'Chile Reef' = hvs[["traits_Chile_reef"]])
Nepean <- list(Global = hvsG[[1]], 'Nepean River' = hvs[["traits_Nepean"]])
Ythan <- list(Global = hvsG[[1]], 'Ythan Estuary' = hvs[["traits_Ythan"]])

compars <- list(Gogo, Miguasha, SantaCruz, BracoMorto, Caribbean, Chile, Nepean, Ythan)

# Define a consistent color palette for the hypervolumes
colors <- c("black", "red") # Generate a unique color for each hypervolume

# Open the PDF device
pdf("/Users/llew0024/Dropbox/Devonian fish/manuscript/ProcB revision/supplementary material/space plots/Sfigure_communityVSglobald.pdf")

# Loop through each comparison in the compars list
for (comparison_name in 1:8) {
  comparison <- compars[[comparison_name]]
  
  # Combine hypervolumes for the current comparison
  combine <- hypervolume_join(comparison)
  names(combine@HVList) <- names(comparison)  # Assign names to the hypervolumes
  
  # Define the pair plot region (top-right half of the page)
  par(fig = c(0, 0.7, 0.3, 1))  # Left, Right, Bottom, Top for the plots
  par(mar = c(4, 4, 2, 1))  # Margins for the plot
  
  # Plot the hypervolumes with the consistent color palette
  plot(combine, 
       show.legend = FALSE, 
       show.data = FALSE,
       main = paste("Comparison:", comparison_name),
       #contour.kde.level = 0.99, # Dynamic title
       show.contour = FALSE,
       col = colors)  # Use the defined color palette for the plots
  
  # Define the legend region (bottom-left part of the page)
  par(fig = c(0, 0.3, 0, 0.3), new = TRUE)  # Left, Right, Bottom, Top for the legend
  par(mar = c(1, 1, 1, 1))  # Minimal margins for the legend
  plot.new()  # Start a new blank plotting area
  
  # Add the legend with the same colors
  legend("center", 
         legend = names(comparison),  # Use the names of the hypervolumes
         col = colors,                # Use the same color palette
         pch = 16,                    # Point symbol
         cex = 0.6)                   # Text size
}

# Close the PDF device
dev.off()

