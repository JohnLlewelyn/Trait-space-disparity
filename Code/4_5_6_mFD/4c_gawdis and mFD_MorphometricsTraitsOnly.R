##Use mFD to apply PCoA and get some functional diversity metrics

#Load packages
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(cowplot)

#Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)


#get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source(paste(filepath,"/code/2_3_tidy data/2_combine_and_tidy_modern_RDS files.R", sep=""), local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

setwd(filepath)

#get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source(paste(filepath,"/code/2_3_tidy data/3_Devonian_tidy.R", sep=""), local = new_env)
dv <- get("mg", envir = new_env)
rm(new_env)

#remove the extra stuff
rm(list= ls()[! (ls() %in% c('fish','dv'))])

#stick it together
fish <- fish[,names(dv)]
fish <- rbind(fish,dv)

#remove Little Rock Lake; too few fish
fish <- fish[fish$community!="traits_Little_Rock_Lake",]

#separate trait data, dropping SL because of inconsistency in how it is measured, and drop mandible because it is highly skewed
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

#combine eye.position front-facing and raised/top of head (only Bothriolepis canadensis has front-facing, but they could also be classified as raised)
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

#cut to morphometric traits)
traits <- traits[,c("logTL","HL","ED", "logPOL", "BD")]
rm(trait, fish, nums)
##get distances using gawdis###################################################################
GD_all <- gawdis(traits) #no issues or warnings 
attr(GD_all,"correls") #contribution are equal
attr(GD_all,"weights") #all have positive weights

qual1 <- quality.fspaces(sp_dist = GD_all, fdendro = "average",maxdim_pcoa = 5,deviation_weighting = c("absolute", "squared"),fdist_scaling = c(TRUE, FALSE)) 

# #check if any species have 0 distance
# duplicates_all <- lapply( traits_cut, function(x){duplicated(x) | duplicated(x, fromLast = TRUE)})
# table(duplicates_all)
# lapply(duplicates_all, table) #no duplicates in terms of trait sets

#check mad index to identify best functional space
# retrieve the functional space associated with minimal quality metric: 
#lapply(qual,function(x){apply(x$quality_fspaces, 2, which.min)}) 
apply(qual1$quality_fspaces,2,which.min) # use all 5
round(qual1$quality_fspaces[1:5, ], 4) #any improvements after 4 or 5 axes are small if any, so use 5 axes

#Alternatively, choose by how much variation explained by each PCoA
eigs <- qual1$details_fspaces$pc_eigenvalues
Esum <- sum(eigs$Eigenvalues)
eigs$perc <- eigs$Eigenvalues/Esum
cbind(rownames(eigs),cumsum(eigs$perc))
plot(rownames(eigs),cumsum(eigs$perc))
print(cumsum(eigs$perc)[4]) #4 eigs  explain 94%
print(cumsum(eigs$perc)[5]) #5 eigs explain 100% - go with that


#position of species
sp_coords1 <- qual1$details_fspaces$sp_pc_coord

#get correlations between traits and axes 
Tcorrs1 <- traits.faxes.cor(sp_tr = traits, sp_faxes_coord = sp_coords1[,paste("PC",1:5, sep="")])
Tcorrs1 <- Tcorrs1[order(Tcorrs1$axis, Tcorrs1$value), ]
#save it
filepath <- "filepath/"
#write.csv(Tcorrs1,paste(filepath,"/tables/PC_trait_correlsHighAccuracy.csv", sep=""), row.names = FALSE)

#calculate functional diversity indices
#need matrix of 1s and 0s where row = community and column = fish species -> the presence_matrix
alpha_FD3 <- alpha.fd.multidim(sp_faxes_coord = sp_coords1[,paste("PC",1:5, sep="")], asb_sp_w = as.matrix(presence_matrix))

#functional diversity according to different metrics
inds3 <- alpha_FD3$functional_diversity_indices  

#and details of distribution
dets3 <- alpha_FD3$details

indsC3 <- inds3[,names(inds3)%in%c("fdis","feve","fric","fmpd",
                                    "fnnd","fdiv","fori","fspe")]


# Convert the data frame from wide to long format
longFM <- function(x) {
  x$site <- row.names(x)
  x$habitat <- ifelse(x$site%in%c("Gogo","traits_Caribbean","traits_Chile_reef"), "reef",
                      ifelse(x$site%in%c("Miguasha","traits_Ythan","Santa_Cruz_Channel"), "estuary","fresh water"))
  x$group <- ifelse(x$site%in%c("Miguasha", "Gogo", "Canowindra"), "Devonian",
                    ifelse(x$site%in%c("Santa_Cruz_Channel","traits_BracoMorto","traits_Caribbean"),"tropical","temperate/subtropical"))
  long_indsC <- pivot_longer(x, 
                             cols = c("fdis","feve","fric","fmpd","fnnd","fdiv","fori","fspe"), 
                             names_to = "metric", 
                             values_to = "value")
}

long_indsC3 <- longFM(indsC3)

#tidy site names
Fnames <- function(x){
  x$site <- gsub("traits_","",x$site)
  x <- as.data.frame(x)
  add1 <- x[x$group%in%"Devonian",] #to make Devonian sites stand out
  x <- rbind(x,add1)
  return(x)}

long_indsC3 <- Fnames(long_indsC3)

#Fix names of metrics
#functional dispersion, functional evenness, functional richness
Mnames <- function(x) {
  x$metric <- ifelse(x$metric=="fdis","dispersion",
                     ifelse(x$metric=="feve", "evenness", 
                            ifelse(x$metric=="fric", "richness",
                                   ifelse(x$metric=="fmpd", "pairwise\ndist.",
                                          ifelse(x$metric=="fnnd", "nearest\nneighbour",
                                                 ifelse(x$metric=="fdiv", "divergence",
                                                        ifelse(x$metric=="fori", "originality","specialization")))))))
  return(data.frame(x))}

long_indsC3 <- Mnames(long_indsC3)

#Create plots in ggplot 
#colour indicates ancient , modern, tropical, temperate
#shape indicates habitat type
#function for plotting the different metrics
plotF <- function(long_indsC1) {
  ggplot(long_indsC1, aes(x = metric, y = value, shape = habitat, color = group)) +
    geom_point(size = 5, alpha = 0.6) +  # Your existing points
    #geom_point(data = long_indsC1[long_indsC1$site == "Little_Rock_Lake", ], aes(x = metric, y = value), color = "black", size = 1, alpha = 0.6, show.legend = FALSE) +  # Points for Little_Rock_Lake
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

metricsPlots3 <- plotF(long_indsC3)


#all traits figure
pdf(paste(filepath,"/figures/metrics_morphometrics.pdf",sep=""),width=18, height=15)
metricsPlots3
dev.off()

#cut to the focal metrics
metrics <- c("richness","nearest\nneighbour","evenness","divergence","specialization")
lc3 <- long_indsC3[long_indsC3$metric%in%metrics,]

#set metrics to factor so can control order
lc3$metric <- factor(lc3$metric, levels = c("richness","nearest\nneighbour","evenness","divergence","specialization"))

Plot3 <- plotF(lc3)

#and plot them
pdf(paste(filepath,"/figures/5metrics_morphometrics.pdf",sep=""),width=12, height=10)
Plot3
dev.off()

saveRDS(Plot3,paste(filepath,"/output/Plot3morphometrics.RDS",sep=""))


