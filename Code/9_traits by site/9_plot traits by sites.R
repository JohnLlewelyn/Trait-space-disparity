# Plot traits by sites
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpattern)

# Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)

# Get modern fish data #together&tidied.rds from combine and tidy RDS files.R
new_env <- new.env()
source("filepath/code/2_3_tidy data/2_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

setwd(filepath)

# Get Devonian fish data # Devonian_traits_tidy.rds from Devonian_tidy.R
new_env <- new.env()
source("filepath/code/2_3_tidy data/3_Devonian_tidy.R", local = new_env)
dv <- get("mg", envir = new_env)
rm(new_env)

# Remove the extra stuff
rm(list= ls()[! (ls() %in% c('fish','dv'))])

# Stick it together
fish <- fish[,names(dv)]
fish <- rbind(fish,dv)

# Remove Little Rock Lake; too few fish
fish <- fish[fish$community!="traits_Little_Rock_Lake",]

# Separate trait data, dropping SL because of inconsistency in how it is measured and mandible because it is highly skewed
trait_names <- c("Species","BodyShapeI","BodyShapeII","TL","HL","ED","POL","BD","PosofMouth","eye.position","spiracle","caudal.fin.shape") #"mandible"
trait <- fish[,trait_names]
trait <- unique(trait)
rownames(trait) <- trait$Species
trait$Species <- NULL
details <- fish[,names(fish)%in%c("Species","community", "habitat", "degsFromEquator")]

# Fix BodyShapeI - fusiform had been split into two groups
trait$BodyShapeI <- as.character(trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("fusi",trait$BodyShapeI), "fusiform",trait$BodyShapeI)
trait$BodyShapeI <- ifelse(grepl("short",trait$BodyShapeI), "shortAndOrDeep",trait$BodyShapeI)

# combine eye.position front-facing and raised/top of head (only Bothriolepis canadensis has front-facing, and they are quite raised)
trait$eye.position <- as.character(trait$eye.position)
trait$eye.position <- ifelse(trait$eye.position=="front-facing","raised/top of head",trait$eye.position)

################################################################################
# Check if the quantitative traits should be log-transformed
nums <- trait[sapply(trait, is.numeric)] 
# Plot histograms
lapply(names(nums), function(col_name) {
  hist(nums[[col_name]], main = col_name)
})
# Make the transformations where needed
# TL
nums$TL <-  log(nums$TL)
names(nums)[names(nums)=="TL"] <- "logTL"
trait$TL <-  log(trait$TL)
names(trait)[names(trait)=="TL"] <- "logTL"
# POL
nums$POL <-  log(nums$POL)
names(nums)[names(nums)=="POL"] <- "logPOL"
trait$POL <-  log(trait$POL)
names(trait)[names(trait)=="POL"] <- "logPOL"
#check distribution again
lapply(names(nums), function(col_name) {
  hist(nums[[col_name]], main = col_name)
})

# Make assemblage matrix
# Add a column for presence (1) for each species in each community
details <- details %>% mutate(presence = 1)

# Transform the data frame into a wide format
presence_matrix <- details %>%
  spread(key = Species, value = presence, fill = 0)

# Remove the community column and community details from the matrix
row.names(presence_matrix) <- presence_matrix$community
cd <- presence_matrix[,1:2]
presence_matrix <- presence_matrix[, !names(presence_matrix)%in%c("community","habitat","degsFromEquator")]

# Need trait detail data frame first, with character columns changed to factor and row names = species
trait_det <- data.frame(trait_name=names(trait) ,trait_type=ifelse(sapply(trait,is.numeric),"Q","N"),trait_weight = 1, fuzzy_name = NA)
traits <- lapply(trait, function(x) if(class(x) == "character") as.factor(x) else x)
traits <- as.data.frame(traits)
rownames(traits) <- rownames(trait)
#check content of the data frames
sp.tr.summary(sp_tr = traits, tr_cat = trait_det)

# Stick traits back with site/details
sp <- traits
traits$Species <- row.names(traits)
BySite <- merge(details[,names(details)%in%c("Species","community","habitat","degsFromEquator")],traits, by = "Species", all.x = TRUE)

# Tidy site names
Fnam <- function(x){
  x$community <- gsub("traits_","",x$community)
  x <- as.data.frame(x)
  return(x)}

BySite <- Fnam(BySite)
BySite$community <- ifelse(BySite$community == "BracoMorto", "Braço Morto\nAcima and Abaixo",
                           ifelse(BySite$community == "Caribbean", "Caribbean Reefs",
                                  ifelse(BySite$community == "Chile_reef", "Chile Reefs", 
                                         ifelse(BySite$community == "Nepean", "Nepean River",
                                                ifelse(BySite$community == "Santa_Cruz_Channel", "Santa Cruz Estuary",
                                                       ifelse(BySite$community == "Ythan","Ythan Estuary", 
                                                              ifelse(BySite$community == "Gogo", "Gogo Reef",
                                                                     ifelse(BySite$community == "Miguasha", "Miguasha Estuary",
                                                                            ifelse(BySite$community == "Canowindra", "Canowindra Billabong", BySite$community)))))))))

# Set metrics to factor so can control order
BySite$community <- factor(BySite$community, levels = c("Gogo Reef", "Miguasha Estuary", "Caribbean Reefs", "Santa Cruz Estuary", "Braço Morto\nAcima and Abaixo", "Canowindra Billabong", "Chile Reefs", "Ythan Estuary", "Nepean River"))

# Add climate zone/age grouping
BySite$degsFromEquator <- ifelse(BySite$degsFromEquator == "?", "Devonian", ifelse(BySite$degsFromEquator %in% c("31.5","34","57"), "temperate/subtropical", "tropical"))
names(BySite)[names(BySite)=="degsFromEquator"] <- "group"

# Adjust groups so colours can match figure 1
BySite$group <- ifelse(BySite$group=="Devonian",BySite$group, BySite$habitat)

# Now make the plots
# Define the original palette
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Convert the palette to RGBA with a specified alpha level (separately for black versus the other colours because black needs to be more transparent)
c1 <- sapply(cbbPalette[1], function(col) {
  rgba <- col2rgb(col) / 255  # Convert hex to normalized RGB
  rgb(rgba[1], rgba[2], rgba[3], alpha = 0.4)  # Re-convert to RGBA with alpha
}, USE.NAMES = FALSE)
c2 <- sapply(cbbPalette[2:4], function(col) {
  rgba <- col2rgb(col) / 255  # Convert hex to normalized RGB
  rgb(rgba[1], rgba[2], rgba[3], alpha = 0.7)  # Re-convert to RGBA with alpha
}, USE.NAMES = FALSE)
cb <- c(c1,c1,c2,c1,c2)

# Names for colours so can match
group_colors <- setNames(cb, c("reef","estuary", "reef", "estuary", "freshwater", "freshwater", "reef", "estuary", "freshwater"))

# which communities get hatching?
hatched_communities <- c("Gogo Reef", "Miguasha Estuary", "Canowindra Billabong")

hatched_communities <- c("Gogo Reef", "Miguasha Estuary", "Canowindra Billabong")

plot_by_site <- function(data, y_var, y_lab, psize,
                         hatched_communities = c("Gogo Reef",
                                                 "Miguasha Estuary",
                                                 "Canowindra Billabong")) {
  # factor levels (so x labels line up as before)
  data$community <- factor(
    data$community,
    levels = c("Gogo Reef","Miguasha Estuary","Caribbean Reefs","Santa Cruz Estuary",
               "Braço Morto\nAcima and Abaixo","Canowindra Billabong",
               "Chile Reefs","Ythan Estuary","Nepean River")
  )
  data$group <- factor(data$habitat)
  data$community_pos <- as.numeric(data$community)
  
  # colours for habitats
  group_colors <- c("reef" = "#E69F00", "estuary" = "#56B4E9", "freshwater" = "#009E73")
  
  # constants for the background split (left = first 5 sites)
  xmin_all <- min(data$community_pos, na.rm = TRUE) - 1
  xmax_all <- max(data$community_pos, na.rm = TRUE) + 1
  left_max <- 5.5
  
  # data for the overlay (Devonian only)
  dat_hat <- subset(data, community %in% hatched_communities)
  
  ggplot(data, aes(x = community_pos, y = .data[[y_var]])) +
    # background halves (use annotate with constants)
    annotate("rect", xmin = xmin_all, xmax = left_max, ymin = -Inf, ymax = Inf,
             fill = rgb((255 + 255)/510, (182 + 255)/510, (193 + 255)/510)) +
    annotate("rect", xmin = left_max, xmax = xmax_all, ymin = -Inf, ymax = Inf,
             fill = rgb(1, 1, 0.5)) +
    
    # 1) points behind
    geom_point(
      aes(color = group),
      position = position_jitter(width = 0.2),
      size = psize, shape = 16, alpha = 1
    ) +
    
    # 2) base violins for all communities (no pattern, semi-transparent fill)
    geom_violin(
      aes(group = community, fill = group),
      alpha = 0.7,
      trim = TRUE,
      scale = "width",   # force consistent width scaling
      colour = NA
    ) +
    
    # 3) overlay stripes ONLY for the three Devonian sites (pattern only, no fill)
    ggpattern::geom_violin_pattern(
      data = dat_hat,
      aes(x = community_pos, y = .data[[y_var]], group = community, fill = group),
      trim = TRUE,
      scale = "width",   # match base
      fill = NA,
      pattern = "stripe",          # or "crosshatch"
      pattern_spacing_unit = "mm",
      pattern_spacing      = 0.02,  # mm distance between lines
      pattern_density      = 0.5,  # ↑ closer to 1 = thicker stripes
      pattern_angle        = 45,
      pattern_fill         = "black",
      pattern_colour       = NA,
      #pattern_alpha        = 0.35,
      colour = NA,
      inherit.aes = FALSE
    ) +
    
    # scales
    scale_fill_manual(values = group_colors, name = "Habitat") +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_x_continuous(breaks = 1:9, labels = levels(data$community)) +
    
    # labels/theme
    labs(y = y_lab, x = "Community") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.box  = "vertical",
      legend.box.just = "left"
    ) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 1))  # solid colour keys
      # no pattern legend added, so you keep just the Habitat legend
    )
}


# Apply to numeric traits
TL <- plot_by_site(BySite, "logTL", "total body length (log-transformed)", 3)
HL <- plot_by_site(BySite, "HL", "head length", 3)
ED <- plot_by_site(BySite, "ED", "eye diameter", 3)
POL <- plot_by_site(BySite, "logPOL", "pre-orbital length (log-transformed)", 3)
BD <- plot_by_site(BySite, "BD", "body depth", 3)
BD2 <- BD

# Make separate ED plot
plot_by_site2 <- function(data, y_var, y_lab, psize) {
  # Ensure 'community' and 'group' are treated as factors
  data$community <- factor(data$community)
  data$group <- factor(data$group)
  # Calculate numeric positions for the communities
  data$community_pos <- as.numeric(data$community)
  # Define color mapping
  group_colors <- c("Devonian" = "#000000", "reef" = "#E69F00", "estuary" = "#56B4E9", "freshwater" = "#009E73") # Define as needed
  # Start plotting
  p <- ggplot(data, aes(x = community_pos, y = .data[[y_var]])) +
    geom_rect(aes(xmin = min(community_pos) - 1, xmax = 5.5, ymin = -Inf, ymax = Inf),
              fill = rgb((255 + 255)/510, (182 + 255)/510, (193 + 255)/510), inherit.aes = FALSE) +
    geom_rect(aes(xmin = 5.5, xmax = max(community_pos) + 1, ymin = -Inf, ymax = Inf),
              fill = rgb(1, 1, 0.5), inherit.aes = FALSE) +
    geom_violin(aes(fill = group, group = community), alpha = 0.7) +
    scale_fill_manual(values = group_colors) +
    geom_point(aes(color = group), # Removed 'shape' from aes()
               position = position_jitter(width = 0.2), size = psize, shape = 16,alpha = 0.5) +  # Specified shape = 16
    scale_color_manual(values = group_colors) +
    theme_classic() +
    labs(y = y_lab, x = "Community", fill = "Habitat/Group", shape = "Habitat") +
    scale_x_continuous(breaks = 1:9, labels = levels(data$community)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2.0)),  # Doubles the size of x-axis text
              axis.text.y = element_text(size = rel(2.0)),  # Doubles the size of y-axis text
              axis.title.x = element_text(size = rel(2.0)),  # Doubles the size of x-axis title
              axis.title.y = element_text(size = rel(2.0)),  # Doubles the size of y-axis title
              legend.title = element_text(size = rel(2.0)),  # Doubles the size of legend titles
              legend.text = element_text(size = rel(2.0)),  # Doubles the size of legend text
              plot.title = element_text(size = rel(2.0)),  # Doubles the size of plot main title
              plot.subtitle = element_text(size = rel(2.0)),axis.ticks.x = element_line(),
              axis.ticks.y = element_line()) +  # Doubles the size of plot subtitle
    guides(fill = guide_legend(override.aes = list(shape = NA)),
           color = "none",
           shape = guide_legend(override.aes = list(color = "black")))
  return(p)
}

ED2 <- plot_by_site2(BySite, "ED", "eye diameter", 3)
ED2 <- ED2 + annotate("text", label = "tropical", x = 1, y = 0.12, angle = 0, size = 8, color = "red") +
       annotate("text", label = "temperate/\nsubtropical", x = 7, y = 0.11, angle = 0, size = 8, color = "black")
  

BD2 <- plot_by_site2(BySite, "BD", "body depth", 3)

# Fix them for plotting
TL <- TL + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + 
  annotate("text", label = "a.", x = 0.8, y = 6.75, size = 16) + annotate("text", label = "tropical", x = 4.5, y = 6.75, angle = 0, size = 12, color = "red", fontface = "bold") +
  annotate("text", label = "temperate/\nsubtropical", x = 7.5, y = 6.6, angle = 0, size = 12, color = "black", fontface = "bold") #axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
HL <- HL + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "b.", x = 0.8, y = 0.46, size = 16)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
ED <- ED + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "c.", x = 0.8, y = 0.12, size = 16)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
POL <- POL + theme(legend.position = "none",axis.text = element_text(size = 30), axis.title.x = element_blank(), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "d.", x = 0.8, y = -1.325, size = 16)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
BD <- BD + theme(legend.position = "none", axis.title.x = element_blank(), axis.text = element_text(size = 30), axis.title = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "e.", x = 0.8, y = 0.64, size = 16)
leg <-BD2 + theme(legend.text = element_text(size = 30),legend.title = element_text(size = 30), legend.key.size = unit(3, "lines"), axis.title.x = element_blank(), axis.text.y = element_text(size = 30), axis.title = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) 
legend_BD2 <- get_legend(leg)

plots_aligned <- cowplot::align_plots(TL, HL, ED, POL, BD, align = 'hv')

# Save the plot as a PNG file
png("filepath/figures/for manuscript/figure 4. numeric_traits_stripes.png", width = 30, height = 40, units = "in", res = 300)
combined_plot <- plot_grid(
  plots_aligned[[1]], plots_aligned[[2]], plots_aligned[[3]],
  plots_aligned[[4]], plots_aligned[[5]], legend_BD2,
  ncol = 2, align = 'hv'
)
print(combined_plot)
dev.off()

##############################
# Now the categorical variables
BySite$BodyShapeI <- as.character(BySite$BodyShapeI)
# Fix one of the categorical variable's levels
BySite$BodyShapeI <- ifelse(BySite$BodyShapeI=="shortAndOrDeep", "short and/or deep", BySite$BodyShapeI)

# define once (or pass in)
hatched_communities <- c("Gogo Reef", "Miguasha Estuary", "Canowindra Billabong")

plot_stacked_bar <- function(data, cat_var) {
  df <- data %>%
    group_by(community, cat = .data[[cat_var]]) %>%
    summarise(n = n(), .groups = "drop")
  
  dev_sites <- hatched_communities
  comm_levels <- levels(factor(df$community))
  # axis_labels <- sapply(comm_levels, function(lbl) {
  #   if (lbl %in% dev_sites) paste0("**<span style='color:black;'>", lbl, "</span>**") else lbl
  # })
  
  ggplot(df, aes(x = community, y = n, fill = cat)) +
    # base bars
    geom_col(position = "fill", colour = NA) +
    
    # striped overlay for Devonian sites only
    ggpattern::geom_col_pattern(
      data = dplyr::filter(df, community %in% dev_sites),
      aes(x = community, y = n, fill = cat),
      pattern = "stripe",           # fixed pattern
      position = "fill",
      colour = NA,
      pattern_spacing_unit = "mm",
      pattern_size_unit    = "mm",
      pattern_spacing      = 0.03,
      pattern_density      = 0.5,
      pattern_angle        = 45,
      pattern_fill         = scales::alpha("black", 1),
      pattern_colour       = NA,
      show.legend = FALSE           # <<< turn off legend contribution
    ) +
    
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Percentage", fill = cat_var) +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 30),
      legend.title = element_blank(),
      legend.key.size = unit(3, "lines"),
      legend.position = "top",
      plot.title = element_text(size = 40, hjust = 0.07, face = "bold"),
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 40),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 40),
      axis.title.y = element_text(size = 40)
    ) +
    #scale_x_discrete(limits = comm_levels, labels = axis_labels) +
    scale_x_discrete(limits = comm_levels, labels = comm_levels) +
    guides(pattern = "none")   # redundant safety — hides pattern legend
}


# Make plots
bs1 <- plot_stacked_bar(BySite, "BodyShapeI") + ggtitle("sagittal body shape  ") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
bs2 <- plot_stacked_bar(BySite, "BodyShapeII") + ggtitle("transverse body shape") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
spir <- plot_stacked_bar(BySite, "spiracle") + ggtitle("spiracle             ") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
tail <- plot_stacked_bar(BySite, "caudal.fin.shape") + ggtitle("caudal fin shape     ") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
mouth <- plot_stacked_bar(BySite, "PosofMouth") + ggtitle("mouth position        ") 
eye <- plot_stacked_bar(BySite, "eye.position") + ggtitle("eye position                   ") 

#####
# # Combine plots with specific row heights
# combined_plots <- plot_grid(
#   plot_grid(bs1, spir, nrow = 1, labels = c("a.", "b."),label_size = 40, label_x = 0.1, label_y = 1.007), # First row with 2 plots
#   plot_grid(bs2, tail, nrow = 1, labels = c("c.", "d."),label_size = 40, label_x = 0.1, label_y = 1.007), # Second row with 2 plots
#   plot_grid(mouth, eye, nrow = 1, labels = c("e.", "f."),label_size = 40, label_x = 0.1, label_y = 1.007), # Third row with 2 plots
#   ncol = 1,
#   rel_heights = c(0.8, 0.83, 1.065) # Adjust these values to set relative row heights
# )

combined_plots <- plot_grid(
  plot_grid(bs1, spir, nrow = 1, labels = c("a.", "b."),
            label_size = 40, label_x = 0.1, label_y = 1.007),
  NULL, # <- spacer row
  plot_grid(bs2, tail, nrow = 1, labels = c("c.", "d."),
            label_size = 40, label_x = 0.1, label_y = 1.007),
  NULL, # <- spacer row
  plot_grid(mouth, eye, nrow = 1, labels = c("e.", "f."),
            label_size = 40, label_x = 0.1, label_y = 1.007),
  ncol = 1,
  rel_heights = c(0.8, 0.05, 0.833, 0.05, 1.084)  # adjust white-space height
)


final_plot <- plot_grid(
  NULL,                # empty space at top
  combined_plots,
  ncol = 1,
  rel_heights = c(0.025, 1)  # 10% of height is white space, adjust as needed
)

# Save figure
png("filepath/figures/for manuscript/figure 5 categorical_traits_stripes.png", width = 30, height = 60, units = "in", res = 300)
#print(combined_plots)
print(final_plot)
dev.off()

# some summary stats
mnsd <- function(x){c(mean=mean(x),standDev=sd(x))}
fish$period <- ifelse(fish$degsFromEquator=="?", "Dev", "mod")
aggregate(fish$ED,by=list(fish$period),FUN=mnsd)
aggregate(fish$BD,by=list(fish$period),FUN=mnsd)

# now % for categorical data
grp<-aggregate(fish$BodyShapeII,by=list(fish$period,fish$BodyShapeII),FUN=length)
grp$perc<-ifelse(grp$Group.1=="Dev", grp$x/80*100,grp$x/436*100)

caud<-aggregate(fish$caudal.fin.shape,by=list(fish$period,fish$caudal.fin.shape),FUN=length)
caud$perc<-ifelse(caud$Group.1=="Dev", caud$x/80*100,caud$x/436*100)

spir<-aggregate(fish$spiracle,by=list(fish$period,fish$spiracle),FUN=length)
spir$perc<-ifelse(spir$Group.1=="Dev", spir$x/80*100,spir$x/436*100)

sag<-aggregate(fish$BodyShapeI,by=list(fish$period,fish$BodyShapeI),FUN=length)
sag$perc<-ifelse(sag$Group.1=="Dev", sag$x/80*100,sag$x/436*100)
