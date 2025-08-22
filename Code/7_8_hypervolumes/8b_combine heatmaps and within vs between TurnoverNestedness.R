# Build heatmaps and within vs between violin plots for turnover and nestedness, and hypervolume plots
library(ggplot2)
library(gridExtra)
library(cowplot)

#Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)

# Get the individual ggplots (generated in Hypervolumes_fixedBandwidths.R)
# Turnover plots
tTime <- readRDS("filepath/output/time_turn_plotJ.rds")
tHabitat <- readRDS("filepath/output/habitat_turn_plotJ.rds")
tClimate <- readRDS("filepath/output/climate_turn_plotJ.rds")
t_plot <- readRDS("filepath/output/turn_plot.rds")

# Create a layout matrix to define the arrangement
layout_matrix <- matrix(
  c(6,7,8,
    1,  2, 5, 
    1,  3, 5,
    1,  4, 5), 
  nrow = 4, ncol = 3, byrow = TRUE)

# Remove some legends and save legend by itself
tTime <- tTime + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
tHabitat <- tHabitat + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
tClimate <- tClimate + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
legend <- get_legend(tTime + theme(legend.position = "right"))

# Make some mostly blank plots for labels
# Create a blank ggplot object
a._plot <- ggplot() + geom_blank() + theme_void() + annotate("text", hjust = 1, x = 100, y = 0, 
  label = "a. heatmap of Jaccard Turnover", size = 10, fontface =2) + 
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10, unit = "pt")) # Increase left margin
b._plot <- ggplot() + geom_blank() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "", size = 15)
empty_plot <- ggplot() + geom_blank() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = " ", size = 15)

combined_turn <- grid.arrange(
  t_plot, tTime, 
  tHabitat, tClimate, legend,a._plot, b._plot, empty_plot,
  nrow = 4,
  ncol = 3,
  layout_matrix = layout_matrix,
  widths = c(4, 1.5, 1), # Adjust the width of the third column as necessary
  heights = c(0.1, 1, 1, 1)
)

ggsave(plot = combined_turn, 
       filename = "filepath/supplementary material/Sfigure # Jaccard Turnover.png",
       height = 20, width =30,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")

# Now plot Jaccard nestedness
# Nestedness plots
nTime <- readRDS("filepath/output/time_nest_plotJ.rds")
nHabitat <- readRDS("filepath/output/habitat_nest_plotJ.rds")
nClimate <- readRDS("filepath/output/climate_nest_plotJ.rds")
nest_plot <- readRDS("filepath/output/nest_plot.rds")

# Remove some legends and save legend by itself
nTime <- nTime + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
nHabitat <- nHabitat + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
nClimate <- nClimate + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
legend <- get_legend(nTime + theme(legend.position = "right"))


a._Jplot <- ggplot() + geom_blank() + theme_void() + annotate("text", hjust = 1, x = 100, y = 0, 
  label = "a. heatmap of Jaccard nestedness", size = 10, fontface =2) + 
  theme(plot.margin = margin(t = 10, r = 110, b = 0, l = 0, unit = "pt")) # Increase left margin


nest_combined <- grid.arrange(
  nest_plot, nTime, 
  nHabitat, nClimate, legend,a._Jplot, b._plot, empty_plot,
  nrow = 4,
  ncol = 3,
  layout_matrix = layout_matrix,
  widths = c(4, 1.5, 1), # Adjust the width of the third column as necessary
  heights = c(0.1, 1, 1, 1)
)


#Save it
ggsave(plot = nest_combined, 
       filename = "filepath/supplementary material/Sfigure 9 Jaccard nestedness.png",
       height = 20, width =30,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")


# compare withing and between
# function for getting mean and sd
mean_sd <- function(x, na.rm = TRUE) {
  c(mean = mean(x, na.rm = na.rm), sd = sd(x, na.rm = na.rm))
}

# groups
dev <- c("Canowindra Billabong","Gogo Reefs", "Miguasha Estuary")
mod <- unique(c(DCentroids$data$Row[!DCentroids$data$Row%in%dev],DCentroids$data$Column[!DCentroids$data$Row%in%dev])) 

#function to select the correct comparisons
btw_rowsF <- function(x) {(x$data$Row %in% dev & x$data$Column %in% mod) |
    (x$data$Row %in% mod & x$data$Column %in% dev)}

# turnover
names(t_plot$data) <- c("Row", "Column", "jaccard")
btwJ <- t_plot$data$jaccard[btw_rowsF(t_plot)]
wiJ <- t_plot$data$jaccard[!btw_rowsF(t_plot)]
lapply(list(btwJ,wiJ),FUN=mean_sd)

# nestedness
names(nest_plot$data) <- c("Row", "Column", "jaccard")
btwJ <- nest_plot$data$jaccard[btw_rowsF(nest_plot)]
wiJ <- nest_plot$data$jaccard[!btw_rowsF(nest_plot)]
lapply(list(btwJ,wiJ),FUN=mean_sd)
