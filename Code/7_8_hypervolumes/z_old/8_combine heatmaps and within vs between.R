#heatmaps combined with within vs between
library(ggplot2)
library(gridExtra)
library(cowplot)

#Adjust file path at"###/"

setwd("###/")

#get the individual ggplots (generated in Hypervolumes_fixedBandwidths.R)
#distance plots
DTime <- readRDS("###/plot resources/dist_time_plot.rds")
DHabitat <- readRDS("###/plot resources/dist_habitat_plot.rds")
DClimate <- readRDS("###/plot resources/dist_climate_plot.rds")
DCentroids <- readRDS("###/plot resources/centroid_dist.rds")

#overlap plots
JTime <- readRDS("###/plot resources/dist_time_plotJ.rds")
JHabitat <- readRDS("###/plot resources/dist_habitat_plotJ.rds")
JClimate <- readRDS("###/plot resources/dist_climate_plotJ.rds")
jac_plot <- readRDS("###/plot resources/jack_plot.rds")

# Create a layout matrix to define the arrangement
layout_matrix <- matrix(
  c(6,7,8,
    1,  2, 5, 
    1,  3, 5,
    1,  4, 5), 
  nrow = 4, ncol = 3, byrow = TRUE)

# remove some legends and save legend by itself
DTime <- DTime + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
DHabitat <- DHabitat + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
DClimate <- DClimate + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
legend <- get_legend(DTime + theme(legend.position = "right"))
#make some mostly blank plots for labels
# Create a blank ggplot object
a._plot <- ggplot() + geom_blank() + theme_void() + annotate("text", hjust = 1, x = 100, y = 0, 
  label = "a. heatmap of centroid distances", size = 10, fontface =2) + 
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10, unit = "pt")) # Increase left margin
b._plot <- ggplot() + geom_blank() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "", size = 15)
empty_plot <- ggplot() + geom_blank() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = " ", size = 15)

combined_dist <- grid.arrange(
  DCentroids, DTime, 
  DHabitat, DClimate, legend,a._plot, b._plot, empty_plot,
  nrow = 4,
  ncol = 3,
  layout_matrix = layout_matrix,
  widths = c(4, 1.5, 1), # Adjust the width of the third column as necessary
  heights = c(0.1, 1, 1, 1)
)

ggsave(plot = combined_dist, 
       filename = "figures/for_ms/figure 3. combine_dist.png",
       height = 20, width =30,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")


#Now the jaccard index
# remove some legends and save legend by itself
JTime <- JTime + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
JHabitat <- JHabitat + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
JClimate <- JClimate + theme(legend.position = "none", plot.margin = margin(t = 10, r = 10, b = 10, l = 60, unit = "pt"))
legend <- get_legend(JTime + theme(legend.position = "right"))


a._Jplot <- ggplot() + geom_blank() + theme_void() + annotate("text", hjust = 1, x = 100, y = 0, 
  label = "a. heatmap of Jaccard index", size = 10, fontface =2) + 
  theme(plot.margin = margin(t = 10, r = 110, b = 0, l = 0, unit = "pt")) # Increase left margin


jacs_combined <- grid.arrange(
  jac_plot, JTime, 
  JHabitat, JClimate, legend,a._Jplot, b._plot, empty_plot,
  nrow = 4,
  ncol = 3,
  layout_matrix = layout_matrix,
  widths = c(4, 1.5, 1), # Adjust the width of the third column as necessary
  heights = c(0.1, 1, 1, 1)
)

ggsave(plot = jacs_combined, 
       filename = "###/figures/for_ms/figure 4. combine_JaccInd.png",
       height = 20, width =30,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")

