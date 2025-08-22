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

setwd("###/")

#get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source("2_3_tidy data/2_Devonian_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

setwd("###/")

#get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source("2_3_tidy data/3_Devonian_tidy_Gogo_Miguasha.R", local = new_env)
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

#Use first 8 PCoAs to build hypervolumes, get the average bandwidths for each axis and then set them as fixed bandwidths (so the same are used for each assemblage)
#PCoA axes are roughly normally distributed, so use Guassian method
set.seed(123) 
start_time <- Sys.time()
hvsprep <- lapply(site_data_frames1, function(x){x <-x[,1:8]})
hvs1 <- lapply(hvsprep, hypervolume_gaussian)
end_time <- Sys.time()
end_time - start_time #1.145 minutes

#compare hypervolumes. Similar for results from mFD?
# Initialize an empty data frame
rdf <- data.frame(site=character(), vol=numeric(), stringsAsFactors=FALSE)
# Loop through each hypervolume in hvs
for(i in 1:length(hvs1)) {
  # Append a new row to rdf
  rdf <- rbind(rdf, data.frame(site=names(hvs1)[i], vol=hvs1[[i]]@Volume, stringsAsFactors=FALSE))
} #dissimilar because different bandwidths used

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

#recalculate hvs with specified bandwidths
bandwidths_fixed <- estimate_bandwidth(hvsprep[[1]], method="fixed", value=mns)
hvs <- hypervolume_gaussian(hvsprep[[1]], kde.bandwidth=bandwidths_fixed)

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

vols <- data.frame(site=character(), vol=numeric(), stringsAsFactors=FALSE)
# Loop through each hypervolume in hvs
for(i in 1:length(hvs)) {
  # Append a new row to rdf
  vols <- rbind(vols, data.frame(site=names(hvsprep)[i], vol=hvs[[i]]@Volume, stringsAsFactors=FALSE))
} #volumes quite similar to mFD richness - time period, habitat, and climate zone patterns are the same

#get distances between all the centroids
names(hvs) <- rownames(presence_matrix)

hypervolume_distances <- combn(names(hvs), 2, function(pairNames) {
  hv1 <- hvs[[pairNames[1]]]
  hv2 <- hvs[[pairNames[2]]]
  distance <- hypervolume_distance(hv1, hv2,type="centroid")
  list(pair = pairNames, distance = distance)
}, simplify = FALSE)

#turn it into a data frame
#Initialize an empty matrix
n <- length(hvs) # Number of hypervolumes
distance_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(distance_matrix) <- names(hvs)
colnames(distance_matrix) <- names(hvs)
#Fill the matrix with distances
for(item in hypervolume_distances) {
  pair <- item$pair
  distance <- item$distance
  distance_matrix[pair[1], pair[2]] <- distance
  distance_matrix[pair[2], pair[1]] <- distance # Assuming distance is symmetric
}
# Convert to data frame
distance_df <- as.data.frame(distance_matrix)

library(reshape2)
library(ggplot2)
# First, add row names as a new column in distance_df
distance_df$RowName <- rownames(distance_df)
# Now, use melt and specify RowName as the id variable
distance_long <- melt(distance_df, id.vars = "RowName", variable.name = "ColumnName", value.name = "Distance")
# Rename columns for clarity
colnames(distance_long) <- c("Row", "Column", "Distance")
distance_long$Column <- as.character(distance_long$Column)
# Fix site names
distance_long$Row <- ifelse(distance_long$Row=="traits_Caribbean", "Caribbean reefs",
                            ifelse(distance_long$Row=="traits_Chile_reef", "Chile reefs",
                                   ifelse(distance_long$Row=="Santa_Cruz_Channel", "Santa Cruz estuary",
                                          ifelse(distance_long$Row=="traits_Ythan","Ythan estuary",
                                                 ifelse(distance_long$Row=="traits_BracoMorto", "Braço Morto\nAcima and Abaixo",
                                                        ifelse(distance_long$Row=="traits_Nepean", "Nepean River", distance_long$Row))))))

distance_long$Column <- ifelse(distance_long$Column=="traits_Caribbean", "Caribbean reefs",
                               ifelse(distance_long$Column=="traits_Chile_reef", "Chile reefs",
                                      ifelse(distance_long$Column=="Santa_Cruz_Channel", "Santa Cruz estuary",
                                             ifelse(distance_long$Column=="traits_Ythan","Ythan estuary",
                                                    ifelse(distance_long$Column=="traits_BracoMorto", "Braço Morto\nAcima and Abaixo",
                                                           ifelse(distance_long$Column=="traits_Nepean", "Nepean River", distance_long$Column))))))

#Fix  names
distance_long$Row[distance_long$Row=="Gogo"] <- "Gogo Reef"
distance_long$Row[distance_long$Row=="Miguasha"] <- "Miguasha Estuary"
distance_long$Row[distance_long$Row=="Santa Cruz estuary"] <- "Santa Cruz Estuary"
distance_long$Row[distance_long$Row=="Ythan estuary"] <- "Ythan Estuary"
distance_long$Row[distance_long$Row=="Chile reefs"] <- "Chile Reefs"
distance_long$Row[distance_long$Row=="Caribbean reefs"] <- "Caribbean Reefs"
distance_long$Row[distance_long$Row=="Braco Morto"] <- "Braço Morto\nAcima and Abaixo"
distance_long$Column[distance_long$Column=="Gogo"] <- "Gogo Reef"
distance_long$Column[distance_long$Column=="Miguasha"] <- "Miguasha Estuary"
distance_long$Column[distance_long$Column=="Santa Cruz estuary"] <- "Santa Cruz Estuary"
distance_long$Column[distance_long$Column=="Ythan estuary"] <- "Ythan Estuary"
distance_long$Column[distance_long$Column=="Chile reefs"] <- "Chile Reefs"
distance_long$Column[distance_long$Column=="Caribbean reefs"] <- "Caribbean Reefs"
distance_long$Column[distance_long$Column=="Braco Morto"] <- "Braço Morto\nAcima and Abaixo"

#fix the order
desired_row_order <- c("Gogo Reef","Miguasha Estuary","Caribbean Reefs","Chile Reefs","Santa Cruz Estuary","Ythan Estuary","Braço Morto\nAcima and Abaixo","Nepean River")  # Replace with your actual row names in the desired order
desired_column_order <- c("Gogo Reef","Miguasha Estuary","Caribbean Reefs","Chile Reefs","Santa Cruz Estuary","Ythan Estuary","Braço Morto\nAcima and Abaixo","Nepean River")  # Replace with your actual column names in the desired order

# Adjust the factor levels in distance_long according to the desired order
distance_long$Row <- factor(distance_long$Row, levels = desired_row_order)
distance_long$Column <- factor(distance_long$Column, levels = desired_column_order)
# Replace NAs
distance_long$Distance <- ifelse(is.na(distance_long$Distance), 0, distance_long$Distance)

# Only plot half
nep <- distance_long[distance_long$Column=="Nepean River",]
bra <- distance_long[distance_long$Column=="Braço Morto\nAcima and Abaixo"&!distance_long$Row%in%c("Nepean River"),]
yth <- distance_long[distance_long$Column=="Ythan Estuary"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River"),]
san <- distance_long[distance_long$Column=="Santa Cruz Estuary"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary"),]
chi <- distance_long[distance_long$Column=="Chile Reefs"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary"),]
car <- distance_long[distance_long$Column=="Caribbean Reefs"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs"),]
mig <- distance_long[distance_long$Column=="Miguasha Estuary"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs","Caribbean Reefs"),]
gog <- distance_long[distance_long$Column=="Gogo Reef"&!distance_long$Row%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs","Caribbean Reefs", "Miguasha Estuary"),]
distance_cut <- rbind(nep,bra,yth,san,chi,car,mig,gog)

# Make distance within versus between violin plots
Dev <- c("Gogo Reef","Miguasha Estuary")
tm <- distance_cut
tm$tm1 <- tm$Row%in%Dev 
tm$tm2 <- tm$Column%in%Dev
tm$WvB <- ifelse(tm$tm1==tm$tm2, "within","between")
tm <- tm[tm$Distance!=0,]
names(tm)[names(tm)=="Distance"] <- "distance"
tm$WvB <- factor(tm$WvB, levels = c("within", "between"))
cb <- c("#009E73","#009E73","#009E73") #c("#009E73","#0072B2","red")
tm$specP <- factor(ifelse(tm$tm1 == TRUE & tm$tm2 == TRUE, "Dev:Dev",ifelse(tm$tm1==FALSE&tm$tm2==FALSE, "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

# Create the plot
# Time period first
dist_time_plot <- ggplot(tm, aes(x = WvB, y = distance, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.2, label = "b. time period", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.2)


saveRDS(dist_time_plot, "plot resources/dist_time_plot.rds")

#same for habitat
hb <- distance_cut
hb$Row<-as.character(hb$Row)
rf <- c("Gogo Reef","Caribbean Reefs","Chile Reefs")
est <- c("Miguasha Estuary", "Santa Cruz Estuary", "Ythan Estuary")
fw <- c("Braço Morto\nAcima and Abaixo","Nepean River")
hb$hb1 <- ifelse(hb$Row%in%rf,"rf",ifelse(hb$Row%in%est, "est", "fw"))
hb$hb2 <- ifelse(hb$Column%in%rf,"rf",ifelse(hb$Column%in%est, "est", "fw"))
hb$Dev <- ifelse(hb$Row%in%Dev&hb$Column%in%Dev,"Dev:Dev",ifelse(!(hb$Row%in%Dev)&!(hb$Column%in%Dev), "mod:mod", "Dev:mod"))
hb<-hb[hb$Distance!=0,]
names(hb)[names(hb)=="Distance"]<-"distance"
hb$WvB <- ifelse(hb$hb1==hb$hb2, "within","between")
hb$WvB <- factor(hb$WvB, levels = c("within", "between"))


cb <- c("#009E73","#009E73","#009E73") #c("#009E73","#0072B2","red")
hb$specP <- factor(ifelse(hb$Dev == "Dev:Dev", "Dev:Dev",ifelse(hb$Dev=="mod:mod", "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

dist_habitat_plot <- ggplot(hb, aes(x = WvB, y = distance, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.2, label = "c. habitat type", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.2)

saveRDS(dist_habitat_plot, "plot resources/dist_habitat_plot.rds")

#now compare climate zones
cl <- distance_cut
cl$Row<-as.character(cl$Row)
trop <- c("Gogo Reef","Caribbean Reefs","Miguasha Estuary", "Santa Cruz Estuary","Braço Morto\nAcima and Abaixo")
cl$cl1 <- ifelse(cl$Row%in%trop,"tropical", "subtropical/temperate")
cl$cl2 <- ifelse(cl$Column%in%trop,"tropical", "subtropical/temperate")
cl$Dev <- ifelse(cl$Row%in%Dev&cl$Column%in%Dev,"Dev:Dev",ifelse(!(cl$Row%in%Dev)&!(cl$Column%in%Dev),"mod:mod", "Dev:mod"))
cl<-cl[cl$Distance!=0,]
names(cl)[names(cl)=="Distance"]<-"distance"
cl$WvB <- ifelse(cl$cl1==cl$cl2, "within","between")
cl$WvB <- factor(cl$WvB, levels = c("within", "between"))
cl$specP <- factor(ifelse(cl$Dev == "Dev:Dev", "Dev:Dev",ifelse(cl$Dev=="mod:mod", "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

dist_climate_plot <- ggplot(cl, aes(x = WvB, y = distance, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.2, label = "d. climate zone", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.2)

saveRDS(dist_climate_plot, "plot resources/dist_climate_plot.rds")

#Create the heatmap
###################
#cut out the bottom (redundant) row
distance_cut <- distance_cut[distance_cut$Distance>0,]
#plot
centroid_dist <- ggplot(distance_cut, aes(x = Row, y = Column, fill = Distance)) +
  geom_tile() + # Creates the tiles for the heatmap
  geom_text(aes(label = sprintf("%.2f", Distance)), color = "black", size = 10) + # Annotates each cell with the distance value
  scale_fill_gradient(low = "white", high = "red", limits = c(NA, NA)) + # Customize the color gradient; adjust as needed
  theme_minimal() + # Minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22), # Rotate x-axis labels for better readability
        axis.text.y = element_text(size = 22),
        axis.title = element_blank(),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.key.size = unit(2, "cm")) + # Remove axis titles for a cleaner look
        labs(fill = "distance\nbetween\ncentroids") # Adjust legend title

saveRDS(centroid_dist, "plot resources/centroid_dist.rds")

pdf("plot resources/centroid_distances_heatmap.pdf", width = 12, height = 11)
centroid_dist
dev.off()

#Overlap
names(hvs) <- c("Gogo","Miguasha","Santa Cruz estuary", "Braço Morto\nAcima and Abaixo", "Caribbean Reefs", "Chile Reefs", "Nepean River", "Ythan estuary")
set.seed(123)
overlap_stats_list <- combn(names(hvs), 2, simplify = FALSE, FUN = function(pairNames) {
  hv1 <- hvs[[pairNames[1]]]
  hv2 <- hvs[[pairNames[2]]]
  # Create the hypervolume set 
  hv_set <- hypervolume_set(hv1, hv2, check.memory = FALSE)
  # Calculate overlap statistics for the set
  set_overlap_stats <- hypervolume_overlap_statistics(hv_set)
  # Naming each item in the list based on the names of the hypervolumes being compared
  list_name <- paste(pairNames[1], pairNames[2], sep = "_vs_")
  # Return a named list for this pair's overlap statistics
  setNames(list(set_overlap_stats), list_name)
})
saveRDS(overlap_stats_list,"plot resources/overlap_stats_list.RDS")

#turn it into a data frame
unList <- function(x){
  x <- c(names(x),unlist(x))
  return(x)
}
overL <- lapply(overlap_stats_list, unList)
overL <- data.frame(do.call("rbind",overL))
names(overL) <- c("comparison","jaccard","sorensen","frac_unique_1","frac_unique_2")
#split the comparison column
split_comparison <- strsplit(overL$comparison, "_vs_")
overL$Place1 <- sapply(split_comparison, `[`, 1)  # Extract first part of each split
overL$Place2 <- sapply(split_comparison, `[`, 2)  # Extract second part of each split
overL$comparison <- NULL
jac <- overL[,names(overL)%in%c("jaccard","Place1","Place2")]
jac$jaccard <- as.numeric(jac$jaccard)
jac2 <- jac[,c("jaccard","Place2","Place1")]
names(jac2) <- c("jaccard","Place1","Place2")
jac <- rbind(jac,jac2)
#Fix  names
jac$Place2[jac$Place2=="Gogo"] <- "Gogo Reef"
jac$Place2[jac$Place2=="Miguasha"] <- "Miguasha Estuary"
jac$Place2[jac$Place2=="Santa Cruz estuary"] <- "Santa Cruz Estuary"
jac$Place2[jac$Place2=="Ythan estuary"] <- "Ythan Estuary"
jac$Place2[jac$Place2=="Chile reefs"] <- "Chile Reefs"
jac$Place2[jac$Place2=="Caribbean reefs"] <- "Caribbean Reefs"
jac$Place1[jac$Place1=="Gogo"] <- "Gogo Reef"
jac$Place1[jac$Place1=="Miguasha"] <- "Miguasha Estuary"
jac$Place1[jac$Place1=="Santa Cruz estuary"] <- "Santa Cruz Estuary"
jac$Place1[jac$Place1=="Ythan estuary"] <- "Ythan Estuary"
jac$Place1[jac$Place1=="Chile reefs"] <- "Chile Reefs"
jac$Place1[jac$Place1=="Caribbean reefs"] <- "Caribbean Reefs"
jac$Place1 <- factor(jac$Place1, levels = desired_column_order)
jac$Place2 <- factor(jac$Place2, levels = desired_row_order)

#jac <- unique(jac)
jaccard_matrix <- xtabs(jaccard ~ Place2 + Place1, data = jac)
jaccard_matrix <- matrix(data = as.vector(jaccard_matrix), 
                         nrow = nrow(jaccard_matrix), 
                         byrow = FALSE, 
                         dimnames = dimnames(jaccard_matrix))

jac_df <- as.data.frame(jaccard_matrix)
jac_long <- melt(jaccard_matrix, id.vars = "RowName", variable.name = "ColumnName", value.name = "jaccard")

saveRDS(jac_long, "plot resources/jac_long.rds")

#only plot half
nep <- jac_long[jac_long$Place2=="Nepean River",]
bra <- jac_long[jac_long$Place2=="Braço Morto\nAcima and Abaixo"&!jac_long$Place1%in%c("Nepean River"),]
yth <- jac_long[jac_long$Place2=="Ythan Estuary"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River"),]
san <- jac_long[jac_long$Place2=="Santa Cruz Estuary"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary"),]
chi <- jac_long[jac_long$Place2=="Chile Reefs"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary"),]
car <- jac_long[jac_long$Place2=="Caribbean Reefs"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs"),]
mig <- jac_long[jac_long$Place2=="Miguasha Estuary"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs","Caribbean Reefs"),]
gog <- jac_long[jac_long$Place2=="Gogo Reef"&!jac_long$Place1%in%c("Braço Morto\nAcima and Abaixo","Nepean River","Ythan Estuary","Santa Cruz Estuary","Chile Reefs","Caribbean Reefs", "Miguasha Estuary"),]
jac_cut <- rbind(nep,bra,yth,san,chi,car,mig,gog)

# Create the plots for overlap within versus between
#time period first
Dev <- c("Gogo Reef","Miguasha Estuary")
tm <- jac_cut
tm$tm1 <- tm$Place2%in%Dev 
tm$tm2 <- tm$Place1%in%Dev
tm$WvB <- ifelse(tm$tm1==tm$tm2, "within","between")
tm <- tm[tm$jaccard!=0,]
names(tm)[names(tm)=="jaccard"] <- "Jaccard_index"
tm$WvB <- factor(tm$WvB, levels = c("within", "between"))
cb <- c("#009E73", "#009E73", "#009E73") #cb <- c("#009E73","#0072B2","red") #"#0072B2"
tm$specP <- factor(ifelse(tm$tm1 == TRUE & tm$tm2 == TRUE, "Dev:Dev",ifelse(tm$tm1==FALSE&tm$tm2==FALSE, "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

dist_time_plotJ <- ggplot(tm, aes(x = WvB, y = Jaccard_index, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.4, label = "b. time period", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.4)

saveRDS(dist_time_plotJ, "plot resources/dist_time_plotJ.rds")

#same for habitat
hb <- jac_cut
rf <- c("Gogo Reef","Caribbean Reefs","Chile Reefs")
est <- c("Miguasha Estuary", "Santa Cruz Estuary", "Ythan Estuary")
fw <- c("Braço Morto\nAcima and Abaixo","Nepean River")
hb$hb1 <- ifelse(as.character(hb$Place2)%in%rf,"rf",ifelse(as.character(hb$Place2)%in%est, "est", "fw"))
hb$hb2 <- ifelse(as.character(hb$Place1)%in%rf,"rf",ifelse(as.character(hb$Place1)%in%est, "est", "fw"))
hb$Dev <- ifelse(hb$Place2%in%Dev&hb$Place1%in%Dev,"Dev:Dev",ifelse(!(hb$Place2%in%Dev)&!(hb$Place1%in%Dev), "mod:mod", "Dev:mod"))
hb<-hb[hb$jaccard!=0,]
names(hb)[names(hb)=="jaccard"]<-"Jaccard_index"
hb$WvB <- ifelse(hb$hb1==hb$hb2, "within","between")
hb$WvB <- factor(hb$WvB, levels = c("within", "between"))
hb$specP <- factor(ifelse(hb$Dev == "Dev:Dev", "Dev:Dev",ifelse(hb$Dev=="mod:mod", "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

dist_habitat_plotJ <- ggplot(hb, aes(x = WvB, y = Jaccard_index, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.4, label = "c. habitat type", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.4)

saveRDS(dist_habitat_plotJ, "plot resources/dist_habitat_plotJ.rds")

#now compare climate zones
cl <- jac_cut
cl$Place2<-as.character(cl$Place2)
trop <- c("Gogo Reef","Caribbean Reefs","Miguasha Estuary", "Santa Cruz Estuary","Braço Morto\nAcima and Abaixo")
cl$cl1 <- ifelse(cl$Place2%in%trop,"tropical", "subtropical/temperate")
cl$cl2 <- ifelse(cl$Place1%in%trop,"tropical", "subtropical/temperate")
cl$Dev <- ifelse(cl$Place2%in%Dev&cl$Place1%in%Dev,"Dev:Dev",ifelse(!(cl$Place2%in%Dev)&!(cl$Place1%in%Dev),"mod:mod", "Dev:mod"))
cl<-cl[cl$jaccard!=0,]
names(cl)[names(cl)=="jaccard"]<-"Jaccard_index"
cl$WvB <- ifelse(cl$cl1==cl$cl2, "within","between")
cl$WvB <- factor(cl$WvB, levels = c("within", "between"))
cl$specP <- factor(ifelse(cl$Dev == "Dev:Dev", "Dev:Dev",ifelse(cl$Dev=="mod:mod", "mod:mod", "Dev:mod")), levels = c("Dev:Dev", "mod:mod", "Dev:mod"))

dist_climate_plotJ <- ggplot(cl, aes(x = WvB, y = Jaccard_index, fill = WvB)) +
  geom_violin(alpha = 0.6) +
  geom_point(aes(color = specP), size = 3, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = cb, guide = FALSE) +  # Set guide to FALSE to remove fill legend
  scale_color_manual(values = c("mod:mod" = "black", "Dev:Dev" = "red", "Dev:mod" = "#F0E442"), 
                     name = "time period\ncombination") +
  theme_classic() +
  labs(y = "Distance between centroids", x = "", fill = "comparison\ntype") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "white", colour = "black", size = 1.5),
        legend.key.size = unit(2, "cm")) +
  guides(color = guide_legend(title = "time period\ncombination")) +  # Ensure color legend is included
  annotate("text", x = min(as.numeric(tm$WvB)) - 0.5, y = 0.4, label = "d. climate zone", hjust = 0, 
           vjust = 0.5, fontface = "bold", size = 10) +
  theme(legend.position = "right") +
  ylim(0, 0.4) 

saveRDS(dist_climate_plotJ, "plot resources/dist_climate_plotJ.rds")

#make heat map
#get rid of distance from themselves
jac_cut <- jac_cut[jac_cut$jaccard>0,]
#plot
jac_plot <- ggplot(jac_cut, aes(x = Place1, y = Place2, fill = jaccard)) +
  geom_tile() + # Creates the tiles for the heatmap
  geom_text(aes(label = sprintf("%.2f", jaccard)), color = "black", size = 10) + # Annotates each cell with the distance value
  scale_fill_gradient(high = "white", low = "red", limits = c(NA, NA)) + # Customize the color gradient; adjust as needed
  theme_minimal() + # Minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22), # Rotate x-axis labels for better readability
        axis.text.y = element_text(size = 22),
        axis.title = element_blank(),
        legend.title = element_text(size = rel(2)),  # Scales legend title
        legend.text = element_text(size = rel(2)),
        legend.key.size = unit(2, "cm"))+
        labs(fill = "Jaccard index") 

saveRDS(jac_plot, "plot resources/jack_plot.rds")

pdf("plot resources/Jaccard_overlap_heatmap_flip.pdf", width = 12, height = 11)
jac_plot
dev.off()

#split into comparisons within versus between same time period so mean and SD can be calculates
mods <- c("Nepean River","Braço Morto\nAcima and Abaixo","Ythan Estuary","Santa Cruz Estuary","Chile Reefs","Caribbean Reefs")
devs <- c("Miguasha Estuary", "Gogo Reef")
jac_mods <- jac_cut[jac_cut$Place2%in%mods&jac_cut$Place1%in%mods,]
jac_devs <- jac_cut[jac_cut$Place2%in%devs&jac_cut$Place1%in%devs,]
jac_diff <- jac_cut[jac_cut$Place2%in%mods&jac_cut$Place1%in%devs,]
jac_same <- rbind(jac_mods,jac_devs)

#function for geng mean and SD
mean_SD<- function(x) {
  z1 <- mean(x)
  z2 <- sd(x)
  return(list(mean=z1,  sd=z2))
}
mean_SD(jac_same$jaccard)
mean_SD(jac_diff$jaccard)

dist_mods <- distance_cut[distance_cut$Column%in%mods&distance_cut$Row%in%mods,]
dist_devs <- distance_cut[distance_cut$Column%in%devs&distance_cut$Row%in%devs,]
dist_diff <- distance_cut[distance_cut$Column%in%mods&distance_cut$Row%in%devs,]
dist_same <- rbind(dist_mods,dist_devs)
mean_SD(dist_same$Distance)
mean_SD(dist_diff$Distance)
