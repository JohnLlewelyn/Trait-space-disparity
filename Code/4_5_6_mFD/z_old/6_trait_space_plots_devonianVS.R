#plot trait spaces

#Adjust file path at"###"

alpha_FD3 <- readRDS("###/Data/data_for_plots/alpha_FD3.RDS")

#functional diversity according to different metrics
inds <- alpha_FD3$functional_diversity_indices #for each community: species richness,  Functional Dispersion, Functional Richness etc

#and details of distribution
dets <- alpha_FD3$details

#plot Miguasha and Gogo compare to contemporary assembalges in FRic metric
indsC <-inds[,names(inds)%in%c("fric","fnnd")]

#scale metrics so they can be plotted together
indsC$fric <-scale(indsC$fric)
indsC$fnnd <-scale(indsC$fnnd)

# Convert the data frame from wide to long format
long_indsC <- lapply(indsC, function(x) {
  x$site <- row.names(x)
  x$habitat <- ifelse(x$site%in%c("Gogo","traits_Caribbean","traits_Chile_reef"), "reef",
                      ifelse(x$site%in%c("Miguasha","traits_Ythan","Santa_Cruz_Channel"), "estuary","fresh water"))
  x$group <- ifelse(x$site%in%c("Miguasha", "Gogo"), "Devonian",
                    ifelse(x$site%in%c("Santa_Cruz_Channel","traits_BracoMorto","traits_Caribbean"),"tropical","temperate/sub-trop"))
  long_indsC <- pivot_longer(x, 
                             cols = c("fnnd", "fric"), 
                             names_to = "metric", 
                             values_to = "value")
  return(long_indsC)})

# Tidy site names
long_indsC <- lapply(long_indsC, function(x){
  x$site <- gsub("traits_","",x$site)
  x <- as.data.frame(x)
  add1 <- x[x$group%in%"Devonian",] #to make Devonian sites stand out
  x <- rbind(x,add1)
  return(x)})

# FIX METRIC NAMES
# Functional dispersion, functional evenness, functional richness
long_indsC <- lapply(long_indsC, function(x) {
  x$metric <- ifelse(x$metric=="fdis","functional dispersion",ifelse(x$metric=="feve", "functional eveness", "functional richness"))
  return(data.frame(x))})

# Create plots in ggplot 
# Colour indicates ancient , modern, tropical, temperate
# Shape indicates habitat type
# Function for plotting the different metrics
plotF <- function(long_indsC1) {
  ggplot(long_indsC1, aes(x = metric, y = value, shape = habitat, color = group, size = 2, alpha = 0.6)) +
    geom_point() +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(size = FALSE, alpha = FALSE) +
    labs(x = "functional diversity metric",
         y = "scaled value")} 

metricPlots <- lapply(long_indsC,plotF)
# Arrange them and plot
figure <- ggarrange(metricPlots[[1]],metricPlots[[2]],metricPlots[[3]],
                    ncol = 2, nrow = 2, heights=c(3,3)) 
pdf("###/figures/metricsV2.pdf",width=14, height=14)
figure
dev.off()

# Plot functional indices; alpha.multidim.plot plots combinations of up to four axes
alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_FD[[1]],
  plot_asb_nm              = c("Miguasha", "Gogo"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
#dev.off()

# Iterate over alpha_FD1 and alpha_FD2
comms <- unique(details$community)
comms <- comms[!comms %in% c("Gogo", "Miguasha")]
comms2 <- gsub("traits_", "", comms)
comms2 <- gsub("_", " ", comms2)
comms2[comms2 == "BracoMorto"] <- "Braco Morto"
all_migVS <- list()
for (alpha_index in 1:length(alpha_FD)) {
  migVS <- list()
  for (i in 1:length(comms)) {
    all_dat <- alpha.multidim.plot(
      output_alpha_fd_multidim = alpha_FD[[alpha_index]],
      plot_asb_nm              = c("Miguasha", comms[i]),
      ind_nm                   = c("fric"),
      
      faxes                    = NULL,
      faxes_nm                 = NULL,
      range_faxes              = c(NA, NA),
      color_bg                 = "grey95",
      shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
      size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
      color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
      shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
      shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
      shape_centroid_fspe      = 23,
      color_centroid_fspe      = "black",
      size_sp_nm               = 3, 
      color_sp_nm              = "black",
      plot_sp_nm               = NULL,
      fontface_sp_nm           = "plain",
      save_file                = FALSE,
      check_input              = TRUE) 
    
    all_dat[1]$fric$PC1_PC2 <- all_dat[1]$fric$PC1_PC2 + ggtitle(comms2[i])
    all_dat[1]$fric$PC1_PC3 <- all_dat[1]$fric$PC1_PC3 + ggtitle(comms2[i])
    all_dat[1]$fric$PC1_PC4 <- all_dat[1]$fric$PC1_PC4 + ggtitle(comms2[i])
    all_dat[1]$fric$PC2_PC3 <- all_dat[1]$fric$PC2_PC3 + ggtitle(comms2[i])
    all_dat[1]$fric$PC3_PC4 <- all_dat[1]$fric$PC3_PC4 + ggtitle(comms2[i])
    
    cc <- list(all_dat[1]$fric$PC1_PC2, all_dat[1]$fric$PC1_PC3, all_dat[1]$fric$PC1_PC4, all_dat[1]$fric$PC2_PC3, all_dat[1]$fric$PC3_PC4)
    migVS[[i]] <- c(cc)
  }
  # Store the migVS list in the all_migVS list
  all_migVS[[alpha_index]] <- migVS
  # Dynamically create filename based on alpha_FD variable
  file_name <- paste0("###/MiguashaVScontemporary_moreData_alpha_FD", alpha_index, ".rds")
  # Save the list to an RDS file
  saveRDS(migVS, file_name)
}


#Manually build legend
datP <- data.frame(Habitat=c("Miguasha","contemporary"), cls=c("#1f968BFF", "#DCE319FF"),vals=c(2,2))
datP <- rbind(datP,datP,datP)
#panel with only legend showing
legendP <- ggplot(datP, aes(x = Habitat, y = vals, fill = cls)) +
  geom_violin(trim = FALSE) +  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white") +
  scale_fill_manual(values = unique(datP$cls),labels = datP$Habitat) +   
  labs(title = "", x = "", y = "", fill="community") +
  theme_minimal()+theme(axis.text=element_blank(),  panel.grid = element_blank()) 
#make a blank panel
blankP <- ggplot(NULL, aes(x = NULL, y = NULL)) +
  geom_blank() +  
  labs(x = NULL, y = NULL) +  
  theme_void() +  
  theme(legend.position = "none")  

#combine the plots in one plot
figureFunction <- function(x) {
  x <- ggarrange( legendP, blankP,blankP,blankP,blankP, 
                  x[[5]][[1]], x[[5]][[2]], x[[5]][[3]], x[[5]][[4]], x[[5]][[5]],
                  x[[2]][[1]], x[[2]][[2]], x[[2]][[3]], x[[2]][[4]], x[[2]][[5]], 
                  x[[6]][[1]], x[[6]][[2]], x[[6]][[3]], x[[6]][[4]], x[[6]][[5]], 
                  x[[1]][[1]], x[[1]][[2]], x[[1]][[3]], x[[1]][[4]], x[[1]][[5]],
                  x[[4]][[1]], x[[4]][[2]], x[[4]][[3]], x[[4]][[4]], x[[4]][[5]],
                  x[[3]][[1]], x[[3]][[2]], x[[3]][[3]], x[[3]][[4]], x[[3]][[5]],
                  ncol = 5, nrow = 7, heights=c(1,3,3,3,3,3,3)) #for labels: labels = c("A","B")
  return(x)
}

figure <- lapply(all_migVS, figureFunction)

pdf("###/PCs Migausha vs Fric.pdf", width = 20, height = 20)
for (i in 1:length(figure)) {
  print(figure[[i]])
}
dev.off()


##same as above but for Gogo#######################################################
#get plots of Gogo Functional Richness versus each of the contemporary communities
# Iterate over alpha_FD1 and alpha_FD2
# Iterate over alpha_FD1 and alpha_FD2
comms <- unique(details$community)
comms <- comms[!comms %in% c("Gogo", "Miguasha")]
comms2 <- gsub("traits_", "", comms)
comms2 <- gsub("_", " ", comms2)
comms2[comms2 == "BracoMorto"] <- "Braco Morto"
all_goVS <- list()
for (alpha_index in 1:length(alpha_FD)) {
  goVS <- list()
  for (i in 1:length(comms)) {
    all_dat <- alpha.multidim.plot(
      output_alpha_fd_multidim = alpha_FD[[alpha_index]],
      plot_asb_nm              = c("Gogo", comms[i]),
      ind_nm                   = c("fric"),
      faxes                    = NULL,
      faxes_nm                 = NULL,
      range_faxes              = c(NA, NA),
      color_bg                 = "grey95",
      shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
      size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
      color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
      alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
      shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
      shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
      shape_centroid_fspe      = 23,
      color_centroid_fspe      = "black",
      size_sp_nm               = 3, 
      color_sp_nm              = "black",
      plot_sp_nm               = NULL,
      fontface_sp_nm           = "plain",
      save_file                = FALSE,
      check_input              = TRUE) 
    all_dat[1]$fric$PC1_PC2 <- all_dat[1]$fric$PC1_PC2 + ggtitle(comms2[i])
    all_dat[1]$fric$PC1_PC3 <- all_dat[1]$fric$PC1_PC3 + ggtitle(comms2[i])
    all_dat[1]$fric$PC1_PC4 <- all_dat[1]$fric$PC1_PC4 + ggtitle(comms2[i])
    all_dat[1]$fric$PC2_PC3 <- all_dat[1]$fric$PC2_PC3 + ggtitle(comms2[i])
    all_dat[1]$fric$PC3_PC4 <- all_dat[1]$fric$PC3_PC4 + ggtitle(comms2[i])
    
    cc <- list(all_dat[1]$fric$PC1_PC2, all_dat[1]$fric$PC1_PC3, all_dat[1]$fric$PC1_PC4, all_dat[1]$fric$PC2_PC3, all_dat[1]$fric$PC3_PC4)
    goVS[[i]] <- c(cc)
  }
  # Store the migVS list in the all_migVS list
  all_goVS[[alpha_index]] <- goVS
  # Dynamically create filename based on alpha_FD variable
  file_name <- paste0("###/GogoVScontemporary_moreData_alpha_FD", alpha_index, ".rds")
  # Save the list to an RDS file
  saveRDS(goVS, file_name)
}


#Manually build legend
datP <- data.frame(Habitat=c("Gogo","contemporary"), cls=c("#1f968BFF", "#DCE319FF"),vals=c(2,2))
datP <- rbind(datP,datP,datP)
#panel with only legend showing
legendP <- ggplot(datP, aes(x = Habitat, y = vals, fill = cls)) +
  geom_violin(trim = FALSE) +  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white") +
  scale_fill_manual(values = unique(datP$cls),labels = datP$Habitat) +   
  labs(title = "", x = "", y = "", fill="community") +
  theme_minimal()+theme(axis.text=element_blank(),  panel.grid = element_blank()) 
#make a blank panel
blankP <- ggplot(NULL, aes(x = NULL, y = NULL)) +
  geom_blank() +  
  labs(x = NULL, y = NULL) +  
  theme_void() +  
  theme(legend.position = "none")  

figureG <- lapply(all_goVS, figureFunction)

pdf("###/PCs Gogo vs Fric.pdf", width = 20, height = 20)
for (i in 1:length(figureG)) {
  print(figureG[[i]])
}
dev.off()
