#plot traits by sites
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ggplot2)

#Adjust file path at"###/"

setwd("###/")

#get modern fish data 
new_env <- new.env()
source("###/2_3_tidy data/2_Devonian_combine_and_tidy_modern_RDS files.R", local = new_env)
fish <- get("traits", envir = new_env)
rm(new_env)

#get Devonian fish data 
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

#stick traits back with site/details
sp <- traits
traits$Species <- row.names(traits)
BySite <- merge(details[,names(details)%in%c("Species","community","habitat","degsFromEquator")],traits, by = "Species", all.x = TRUE)

#tidy site names
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
                                                                     ifelse(BySite$community == "Miguasha", "Miguasha Estuary", BySite$community))))))))

#set metrics to factor so can control order
BySite$community <- factor(BySite$community, levels = c("Gogo Reef", "Miguasha Estuary", "Caribbean Reefs", "Santa Cruz Estuary", "Braço Morto\nAcima and Abaixo", "Chile Reefs", "Ythan Estuary", "Nepean River"))

#add climate zone/age grouping
BySite$degsFromEquator <- ifelse(BySite$degsFromEquator == "?", "Devonian", ifelse(BySite$degsFromEquator %in% c("31.5","34","57"), "temperate/subtropical", "tropical"))
names(BySite)[names(BySite)=="degsFromEquator"] <- "group"

#adjust groups so colours can match figure 1
BySite$group <- ifelse(BySite$group=="Devonian",BySite$group, BySite$habitat)

#Now make the plots
# Define the original palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Convert the palette to RGBA with a specified alpha level (separately for black versus the other colours because black needs to be more transparent)
c1 <- sapply(cbbPalette[1], function(col) {
  rgba <- col2rgb(col) / 255  # Convert hex to normalized RGB
  rgb(rgba[1], rgba[2], rgba[3], alpha = 0.4)  # Re-convert to RGBA with alpha
}, USE.NAMES = FALSE)
c2 <- sapply(cbbPalette[2:4], function(col) {
  rgba <- col2rgb(col) / 255  # Convert hex to normalized RGB
  rgb(rgba[1], rgba[2], rgba[3], alpha = 0.7)  # Re-convert to RGBA with alpha
}, USE.NAMES = FALSE)
cb <- c(c1,c1,c2,c2)

#names for colours so can match
group_colors <- setNames(cb, c("Devonian","Devonian", "reef", "estuary", "freshwater","reef", "estuary", "freshwater"))

plot_by_site <- function(data, y_var, y_lab, psize) {
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
    geom_violin(aes(fill = group, group = community), alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    geom_point(aes(color = group), # Removed 'shape' from aes()
               position = position_jitter(width = 0.2), size = psize, shape = 16,alpha = 0.5) +  # Specified shape = 16
    scale_color_manual(values = group_colors) +
    theme_classic() +
    labs(y = y_lab, x = "Community", fill = "Habitat/Group", shape = "Habitat") +
    scale_x_continuous(breaks = 1:8, labels = levels(data$community)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(override.aes = list(shape = NA)),
           color = "none",
           shape = guide_legend(override.aes = list(color = "black")))
  return(p)
}

#apply to numeric traits
TL <- plot_by_site(BySite, "logTL", "total body length (log-transformed)", 6)
HL <- plot_by_site(BySite, "HL", "head length", 6)
ED <- plot_by_site(BySite, "ED", "eye diameter", 6)
POL <- plot_by_site(BySite, "logPOL", "pre-orbital length (log-transformed)", 6)
BD <- plot_by_site(BySite, "BD", "body depth", 6)
BD2 <- BD

#make separate ED plot
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
    geom_violin(aes(fill = group, group = community), alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    geom_point(aes(color = group), # Removed 'shape' from aes()
               position = position_jitter(width = 0.2), size = psize, shape = 16,alpha = 0.5) +  # Specified shape = 16
    scale_color_manual(values = group_colors) +
    theme_classic() +
    labs(y = y_lab, x = "Community", fill = "Habitat/Group", shape = "Habitat") +
    scale_x_continuous(breaks = 1:8, labels = levels(data$community)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2.0)),  # Doubles the size of x-axis text
              axis.text.y = element_text(size = rel(2.0)),  # Doubles the size of y-axis text
              axis.title.x = element_text(size = rel(2.0)),  # Doubles the size of x-axis title
              axis.title.y = element_text(size = rel(2.0)),  # Doubles the size of y-axis title
              legend.title = element_text(size = rel(2.0)),  # Doubles the size of legend titles
              legend.text = element_text(size = rel(2.0)),  # Doubles the size of legend text
              plot.title = element_text(size = rel(2.0)),  # Doubles the size of plot main title
              plot.subtitle = element_text(size = rel(2.0))) +  # Doubles the size of plot subtitle
    guides(fill = guide_legend(override.aes = list(shape = NA)),
           color = "none",
           shape = guide_legend(override.aes = list(color = "black")))
  return(p)
}

ED2 <- plot_by_site2(BySite, "ED", "eye diameter", 4)
ED2 <- ED2 + annotate("text", label = "tropical", x = 1, y = 0.12, angle = 0, size = 8, color = "red") +
       annotate("text", label = "temperate/\nsubtropical", x = 7, y = 0.11, angle = 0, size = 8, color = "black")
  

BD2 <- plot_by_site2(BySite, "BD", "body depth", 4)

#fix them for plotting
TL <- TL + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + 
  annotate("text", label = "a.", x = 0.8, y = 6.75, size = 14) + annotate("text", label = "tropical", x = 4.5, y = 6.75, angle = 0, size = 12, color = "red", fontface = "bold") +
  annotate("text", label = "temperate/\nsubtropical", x = 7.5, y = 6.6, angle = 0, size = 12, color = "black", fontface = "bold") #axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
HL <- HL + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "b.", x = 0.8, y = 0.41, size = 14)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
ED <- ED + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "c.", x = 0.8, y = 0.12, size = 14)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
POL <- POL + theme(legend.position = "none",axis.text = element_text(size = 30), axis.title.x = element_blank(), axis.title.y = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "d.", x = 0.8, y = -1.325, size = 14)#axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + annotate("text", label = "a)", x = 0.6, y = 2.5, size = 7)
BD <- BD + theme(legend.position = "none", axis.title.x = element_blank(), axis.text = element_text(size = 30), axis.title = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) + annotate("text", label = "e.", x = 0.8, y = 0.64, size = 14)
leg <-BD2 + theme(legend.text = element_text(size = 30),legend.title = element_text(size = 30), legend.key.size = unit(3, "lines"), axis.title.x = element_blank(), axis.text.y = element_text(size = 30), axis.title = element_text(size = 30),plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
legend_BD2 <- get_legend(leg)

#plots_aligned <- cowplot::align_plots(TL, HL, ED, POL, BD, legend_BD2, align = 'hv')
plots_aligned <- cowplot::align_plots(TL, HL, ED, POL, BD, align = 'hv')

# #and plot them
# pdf("figures/for_ms/figure 5 numeric_traits.pdf", width = 30, height = 40)
# combined_plot <- plot_grid(
#   plots_aligned[[1]], plots_aligned[[2]], plots_aligned[[3]],
#   plots_aligned[[4]], plots_aligned[[5]], legend_BD2,
#   ncol = 2, align = 'hv'
# )
# print(combined_plot)
# dev.off()

# Save the plot as a PNG file
png("###/figure_5_numeric_traits.png", width = 30, height = 40, units = "in", res = 300)
combined_plot <- plot_grid(
  plots_aligned[[1]], plots_aligned[[2]], plots_aligned[[3]],
  plots_aligned[[4]], plots_aligned[[5]], legend_BD2,
  ncol = 2, align = 'hv'
)
print(combined_plot)
dev.off()

##############################
#now the categorical variables
BySite$BodyShapeI <- as.character(BySite$BodyShapeI)
#fix one of the categorical variable's levels
BySite$BodyShapeI <- ifelse(BySite$BodyShapeI=="shortAndOrDeep", "short and/or deep", BySite$BodyShapeI)

#plot function for categorical variables
plot_stacked_bar <- function(data, cat_var) {
  # Calculate the proportions
  data <- data %>%
    group_by(community, .data[[cat_var]]) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(total = sum(n), percentage = n / total * 100) %>%
    ungroup()

  # Plot
  ggplot(data, aes(x = community, y = percentage, fill = .data[[cat_var]])) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +
    labs(y = "Percentage", fill = cat_var) +
    theme(legend.text = element_text(size = 30),legend.title = element_blank(), legend.key.size = unit(3, "lines"),
          legend.margin = margin(1,1,1,1),
          legend.box.margin = margin(10,1,1,1),
          legend.position="top",
          plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 40),axis.title.x = element_blank(),
          axis.text.y = element_text(size = 40), axis.title.y = element_text(size = 40))
}

# make plots
bs1 <- plot_stacked_bar(BySite, "BodyShapeI") + ggtitle("sagittal body shape") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
bs2 <- plot_stacked_bar(BySite, "BodyShapeII") + ggtitle("transverse body shape") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
spir <- plot_stacked_bar(BySite, "spiracle") + ggtitle("spiracle") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
tail <- plot_stacked_bar(BySite, "caudal.fin.shape") + ggtitle("caudal fin shape") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())
mouth <- plot_stacked_bar(BySite, "PosofMouth") + ggtitle("mouth position") 
eye <- plot_stacked_bar(BySite, "eye.position") + ggtitle("eye position") 

#####
plot_stacked_bar <- function(data, cat_var) {
  # Your existing code up to theme()...
  theme(legend.text = element_text(size = 20),legend.title = element_blank(), legend.key.size = unit(3, "lines"),
        plot.title = element_text(size = 30, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),axis.title.x = element_blank(),axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
    guides(fill = "none") # Add this to remove the legend
}


# Combine plots with specific row heights
combined_plots <- plot_grid(
  plot_grid(bs1, spir, nrow = 1, labels = c("a.", "b."),label_size = 40, label_x = 0.1), # First row with 2 plots
  plot_grid(bs2, tail, nrow = 1, labels = c("c.", "d."),label_size = 40, label_x = 0.1), # Second row with 2 plots
  plot_grid(mouth, eye, nrow = 1, labels = c("e.", "f."),label_size = 40, label_x = 0.1), # Third row with 2 plots
  ncol = 1,
  rel_heights = c(0.8, 0.84, 1.0) # Adjust these values to set relative row heights
)


#save figure
png("###/figure 6 categorical_traits.png", width = 30, height = 50, units = "in", res = 300)
print(combined_plots)
dev.off()


