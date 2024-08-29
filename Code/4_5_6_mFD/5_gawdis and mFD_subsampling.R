#Use mFD to apply PCoA and get some functional richness metrics, controlling for species diversity by sub-sampling populations
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(grid)
library(gridExtra)
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

# calculate functional diversity indices
# need matrix of 1s and 0s where row = community and column = fish species -> the presence_matrix
# loop to sub-sample 10 species from each assemblage and save as a list
set.seed(89)
res <- list()
for(i in 1:100) {
  pm <- presence_matrix[row.names(presence_matrix) != "traits_Little_Rock_Lake", ]
  for (j in 1:nrow(pm)) {
    ones_indices <- which(pm[j, ] == 1)
    if (length(ones_indices) > 10) {
      shuffled_indices <- sample(ones_indices)
      switch_indices <- shuffled_indices[11:length(shuffled_indices)]
      pm[j, switch_indices] <- 0
    }
  }
  #cols to keep
  kp <- which(colSums(pm)>0)
  pm <- pm[,kp]
  res[[i]] <- pm
}

#now get metrics for these subs-samples
appA <- function(x) { alpha.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:8, sep="")], asb_sp_w = as.matrix(x))}

sct <- list()
for (i in 1:3){
  sc <- sp_coords[[i]]
  sct[[i]] <- lapply(res,function(x) { alpha.fd.multidim(sp_faxes_coord = sc[,paste("PC",1:7, sep="")], asb_sp_w = as.matrix(x))})
  }
saveRDS(sct, "###/sct_sub_sampled metrics.RDS")

start_time <- Sys.time()
alpha_FD1 <- lapply(res, appA)
end_time <- Sys.time()
end_time - start_time #takes a bit of time
saveRDS(alpha_FD1, "###/sub_sampled metrics.RDS")

#functional diversity according to different metrics
indsc <- lapply(sct, function(outer_list) {
  lapply(outer_list, function(inner_list) {
    inner_list$functional_diversity_indices
  })
})

inds1 <- lapply(alpha_FD1, function(x) x$functional_diversity_indices)

#and details of distribution
dets <- lapply(sct, function(outer_list) {
  lapply(outer_list, function(inner_list) {
    inner_list$details
  })
})
dets1 <- lapply(alpha_FD1, function(x) x$details)

indsct <- lapply(indsc, function(outer_list) {
  lapply(outer_list, function(x){x<-x[,names(x)%in%c("fdis","feve","fric","fmpd",
                                                     "fnnd","fdiv","fori","fspe")]})
})
  
indsC <- lapply(inds1, function(x){x<-x[,names(x)%in%c("fdis","feve","fric","fmpd",
                                                      "fnnd","fdiv","fori","fspe")]})

#now get the means
indsct <-  lapply(indsct, function(outer_list) {
  lapply(outer_list, function(x){sites <- rownames(x)
                                x<-cbind(x,sites)
                                })
})

indsct <- lapply(indsct, function(list_of_dfs) {
    combined_df <- do.call(rbind, list_of_dfs)})

indsct <- lapply(indsct, function(x){x$sites <- gsub('[0-9]+', '', x$sites)
return(x)})

indsct <- lapply(indsct, function(x){aggregate(x[,1:8], by = list(x$sites), mean)})

indsC <- do.call("rbind",indsC)
indsC$sites <- row.names(indsC)
indsC$sites <- gsub('[0-9]+', '', indsC$sites)
indsC <- aggregate(indsC[,1:8], by = list(indsC$sites), mean)

indsC1 <- indsC

#fix names
indsct <- lapply(indsct, function(x){
  names(x)<- c("site","fdis","fmpd","fnnd","feve","fric","fdiv","fori","fspe")
  return(x)
})
names(indsC1)[names(indsC1)=="Group.1"] <- "site"

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

long_indsct <- lapply(indsct, longFM)
long_indsC1 <- longFM(indsC1)

#tidy site names
Fnames <- function(x){
  x$site <- gsub("traits_","",x$site)
  x <- as.data.frame(x)
  add1 <- x[x$group%in%"Devonian",] #to make Devonian sites stand out
  x <- rbind(x,add1)
  return(x)}

long_indsct <- lapply(long_indsct,Fnames)
long_indsC1 <- Fnames(long_indsC1)

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

long_indsct <- lapply(long_indsct, Mnames)
long_indsC1 <- Mnames(long_indsC1)

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

combo_plots <- lapply(long_indsct, plotF)
metricsPlots1 <- plotF(long_indsC1)

#all traits figure
pdf("###/metrics_all_traits_controlDiversity.pdf",width=12, height=10)
metricsPlots1
dev.off()

#in case it crashes, because it takes a while to get here
save(list = ls(all.names = TRUE), file = "###/my_workspace_compNOunbalanced_subsampling.RData")

#cut to the three focal metrics
metrics <- c("richness","nearest\nneighbour","specialization")
lc1 <- long_indsC1[long_indsC1$metric%in%metrics,]
lc1$metric <- factor(lc1$metric, levels = c("richness", "nearest\nneighbour", "specialization"))
Plot1 <- plotF(lc1)

lc <- lapply(long_indsct,function(x){x <- x[x$metric%in%metrics,]
              x$metric <- factor(x$metric, levels = c("richness", "nearest\nneighbour", "specialization"))
  return(x)})
plotCombo <- lapply(lc,plotF)

pdf("###/3metrics_all_traits_controlDiversity.pdf",width=12, height=10)
Plot1
dev.off()

saveRDS(Plot1,"###/Plot1.RDS")

#combine the combos in a single plot
pdf("###/figure S3. 3metrics_traitCombos_controlDiversity.pdf",width=18, height=10)
final_plot <- plot_grid(
  plotCombo[[1]] + theme(legend.position = "none"),
  plotCombo[[2]] + theme(legend.position = "none"),
  plotCombo[[3]],
  ncol = 3, nrow = 1, rel_widths = c(3, 3, 5)  # Specify the relative widths
)
print(final_plot)
dev.off()
