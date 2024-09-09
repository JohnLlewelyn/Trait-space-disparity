##Use mFD to apply PCoA and get some functional richness metrics, including all species diversity

library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(cowplot)

#Adjust file path at"###"

setwd("###")

#get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source("###/2_3_tidy data/2_Devonian_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

#get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source("###/2_3_tidy data/3_Devonian_tidy_Gogo_Miguasha.R", local = new_env)
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

 #figure out which traits to keep - must have no negative weights (but unbalanced distribution is okay - indicates distribution of values that is heavily skewed or that these traits have many identical values for most of the species#
 num_traits <- ncol(traits)
 max_combinations <- vector("list", num_traits)
 max_no_warning <- 0
 max_combination_names <- list()

# Function to apply gawdis and catch warnings
safe_gawdis <- function(combination) {
  warnings <- NULL
  result <- withCallingHandlers(
    gawdis(traits[, combination, drop = FALSE]),
    warning = function(w) warnings <<- c(warnings, w$message)
  )
  list(result = result, warnings = warnings)
}

# Initialize max_combinations for each possible length
for (i in 8:ncol(traits)) {
  max_combinations[[i]] <- list()
}

# Loop over all combinations starting from 8 - it takes a while
for (i in 8:ncol(traits)) {
  combinations <- combn(num_traits, i, simplify = FALSE)

  for (combination in combinations) {
    result <- safe_gawdis(combination)

    if (is.null(result$warnings)) {  # No warnings
      combination_length <- length(combination)
      if (combination_length > max_no_warning) {
        max_no_warning <- combination_length
        max_combination_names <- list(colnames(trait)[combination])
      } else if (combination_length == max_no_warning) {
        max_combination_names <- c(max_combination_names, list(colnames(trait)[combination]))
      }
    }
  }
}

# max_combination_names contains the names of traits in the models with the most traits that didn't produce warnings
max_combination_names

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
qual <- lapply(gd, quality.fspaces,fdendro = "average",maxdim_pcoa = 10,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE))
qual1 <- quality.fspaces(sp_dist = GD_all, fdendro = "average",maxdim_pcoa = 10,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 

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
saveRDS(species_coords, "###/data/PCoAs/species_coordinates.RDS")

#see correlations between traits and axes #can handle up to 10 traits, produces data frames and plots
Tcorrs <- mapply(function(x,y){
  obs<-traits.faxes.cor(sp_tr = x, sp_faxes_coord = y[,paste("PC",1:7, sep="")], plot = TRUE)}, 
  x=traits_cut, y=sp_coords)
           
#get correlations between traits and axes 
Tcorrs1 <- traits.faxes.cor(sp_tr = traits, sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")])
Tcorrs1 <- Tcorrs1[order(Tcorrs1$axis, Tcorrs1$value), ]
#save it
write.csv(Tcorrs1,"###/PC_trait_correls.csv", row.names = FALSE)

#calculate functional diversity indices
#need matrix of 1s and 0s where row = community and column = fish species -> the presence_matrix
alpha_FD <- lapply(sp_coords, function(x){alpha.fd.multidim(sp_faxes_coord = x[,paste("PC",1:7, sep="")], asb_sp_w = as.matrix(presence_matrix))}) #can scale (scaling = TRUE) so values are between 0 and 1, but different indices are squashed into different portions of this range; can instead set mean to 0 and sd to 1
alpha_FD3 <- alpha.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")], asb_sp_w = as.matrix(presence_matrix))

#functional diversity according to different metrics
inds <- lapply(alpha_FD, function(x){x$functional_diversity_indices}) #for each community: species richness,  Functional Dispersion, Functional Richness etc
inds3 <- alpha_FD3$functional_diversity_indices  

#and details of distribution
dets <- lapply(alpha_FD, function(x){x$details})
dets3 <- alpha_FD3$details

#subset to metrics of interest
indsC <- lapply(inds, function(x){x<-x[,names(x)%in%c("fdis","feve","fric","fmpd",
                                                      "fnnd","fdiv","fori","fspe")]})

indsC3 <- inds3[,names(inds3)%in%c("fdis","feve","fric","fmpd",
                                    "fnnd","fdiv","fori","fspe")]


# Convert the data frame from wide to long format
longFM <- function(x) {
  x$site <- row.names(x)
  x$habitat <- ifelse(x$site%in%c("Gogo","traits_Caribbean","traits_Chile_reef"), "reef",
                      ifelse(x$site%in%c("Miguasha","traits_Ythan","Santa_Cruz_Channel"), "estuary","fresh water"))
  x$group <- ifelse(x$site%in%c("Miguasha", "Gogo"), "Devonian",
                    ifelse(x$site%in%c("Santa_Cruz_Channel","traits_BracoMorto","traits_Caribbean"),"tropical","temperate/subtropical"))
  long_indsC <- pivot_longer(x, 
                             cols = c("fdis","feve","fric","fmpd","fnnd","fdiv","fori","fspe"), 
                             names_to = "metric", 
                             values_to = "value")
}

long_indsC <- lapply(indsC, longFM)
long_indsC3 <- longFM(indsC3)

#tidy site names
Fnames <- function(x){
  x$site <- gsub("traits_","",x$site)
  x <- as.data.frame(x)
  add1 <- x[x$group%in%"Devonian",] #to make Devonian sites stand out
  x <- rbind(x,add1)
  return(x)}

long_indsC <- lapply(long_indsC, Fnames)
long_indsC3 <- Fnames(long_indsC3)

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

long_indsC <- lapply(long_indsC, Mnames)
long_indsC3 <- Mnames(long_indsC3)

#Create plots in ggplot 
#colour indicates ancient , modern, tropical, temperate
#shape indicates habitat type
#function for plotting the different metrics
plotF <- function(long_indsC1) {
  ggplot(long_indsC1, aes(x = metric, y = value, shape = habitat, color = group)) +
    geom_point(size = 5, alpha = 0.6) +  # Your existing points
    geom_point(data = long_indsC1[long_indsC1$site == "Little_Rock_Lake", ], aes(x = metric, y = value), color = "black", size = 1, alpha = 0.6, show.legend = FALSE) +  # Points for Little_Rock_Lake
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

metricPlots <- lapply(long_indsC,plotF)
metricsPlots3 <- plotF(long_indsC3)
#METRICS SHOWING CONSISTENT PATTERNS: nearest neighbour, originality, richness, and specialization
#cut the combo traits results to the 3 consistent results (dropping originality because it is similar to specialisation)
long_indsCc <- lapply(long_indsC, function(x) x <- x[x$metric%in%c("richness","nearest\nneighbour","specialization"),])
metricPlots4 <- lapply(long_indsCc,plotF)

#remove legends from some of them
metricPlots[[1]] <- metricPlots[[1]] + theme(legend.position = "none")
metricPlots[[2]] <- metricPlots[[2]] + theme(legend.position = "none")
metricPlots[[3]] <- metricPlots[[3]] + theme(legend.position = "none")

#arrange them and plot
figure <- ggarrange(metricPlots[[1]],metricPlots[[2]],metricPlots[[3]],
ncol = 3, nrow = 1, heights=c(3), widths = c(3,3,3.3)) 
pdf("###/figure S2. 3combos_balanced_allspecies_alltraits.pdf",width=20, height=10)
figure
dev.off()

#all traits figure
pdf("###/metrics_all_traitsEig5.pdf",width=12, height=10)
metricsPlots3
dev.off()

#cut to the focal metrics
metrics <- c("richness","originality","specialization","nearest\nneighbour")
lc3 <- long_indsC3[long_indsC3$metric%in%metrics,]

#set metrics to factor so can control order
lc3$metric <- factor(lc3$metric, levels = c("richness", "nearest\nneighbour", "originality", "specialization"))

Plot3 <- plotF(lc3)

#and plot them
pdf("###/4metrics_all_traits.pdf",width=12, height=10)
Plot3
dev.off()

saveRDS(Plot3,"###/data/data_for_plots/Plot3.RDS")

#and with the 3 combo of traits
long_indsCc <- lapply(long_indsCc, function(x) { x$metric<-
  factor(x$metric, levels = c("richness", "nearest\nneighbour", "originality", "specialization"))
         return(x)})
metricPlots4 <- lapply(long_indsCc,plotF)
#remove legends from some of them
metricPlots4[[1]] <- metricPlots4[[1]] + theme(legend.position = "none")
metricPlots4[[2]] <- metricPlots4[[2]] + theme(legend.position = "none")

figure <- ggarrange(metricPlots4[[1]],metricPlots4[[2]],metricPlots4[[3]],
                    ncol = 3, nrow = 1, heights=c(3), widths = c(3, 3, 5)) 
pdf("###/figure S2. 3combos_balanced_allspecies_alltraits.pdf",width=18, height=10)
figure
dev.off()

###do it also with sub-sampling, to control for species diversity - see script gawdis and mFD_subsampling.R### ##################
##
##
###get the centroids (no need to sub-ample for this)################################################################################
cr <- function(x, n_coords) {
  coords <- x$asb_G_coord
  # Dynamically bind the specified number of coords
  cents_list <- lapply(1:n_coords, function(i) coords[[i]])
  cents <- do.call(rbind, cents_list)
  cents <- data.frame(cents)
  cents$site <- names(coords)[1:n_coords]
  cents$site <- gsub("traits_", "", cents$site)
  cents$habitat <- ifelse(cents$site %in% c("Gogo", "Caribbean", "Chile_reef"), "reef",
                          ifelse(cents$site %in% c("Miguasha", "Ythan", "Santa_Cruz_Channel"), "estuary", "fresh water"))
  cents$group <- ifelse(cents$site %in% c("Miguasha", "Gogo"), "Devonian",
                        ifelse(cents$site %in% c("Santa_Cruz_Channel", "BracoMorto", "Caribbean"), "tropical", "temperate/subtropical"))
  return(cents)
}

centroids <- lapply(dets, cr, n_coords = 7)
centroids3 <- cr(dets3,8)

#save as table
cent3_round <- centroids3
cent3_round <- cent3_round %>%
  mutate_if(is.numeric, round, digits = 3)
write.csv(cent3_round,"###/PC_trait_centroids.csv", row.names = FALSE)

