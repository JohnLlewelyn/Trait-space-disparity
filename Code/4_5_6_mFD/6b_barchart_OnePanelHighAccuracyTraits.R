#make barchart of mFD results
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(grid)
library(gridExtra)
library(patchwork)
library(colorBlindness)

#Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)

Plot3 <- readRDS("filepath/output/Plot3HighAccuracyTraits.RDS")
d3 <- Plot3$data

#assign data 
df <- d3 
df$metric <- as.character(df$metric)

df$metric <- ifelse(df$metric=="richness", "functional richness", 
                    ifelse(df$metric=="nearest\nneighbour", "nearest neighbour",df$metric))
                          

df$site <- ifelse(df$site == "BracoMorto", "Braço Morto\nAcima and Abaixo",
                           ifelse(df$site == "Caribbean", "Caribbean Reefs",
                                  ifelse(df$site == "Chile_reef", "Chile Reefs", 
                                         ifelse(df$site == "Nepean", "Nepean River",
                                                ifelse(df$site == "Santa_Cruz_Channel", "Santa Cruz Estuary",
                                                       ifelse(df$site == "Ythan","Ythan Estuary", df$site))))))

df$site[df$site=="Gogo"] <- "Gogo Reef"
df$site[df$site=="Miguasha"] <- "Miguasha Estuary"
df$site[df$site=="Canowindra"] <- "Canowindra Billabong"


# Define the order of sites within each group
site_order <- c("Gogo Reef","Miguasha Estuary", "Caribbean Reefs", "Santa Cruz Estuary", "Braço Morto\nAcima and Abaixo","","Canowindra Billabong","Chile Reefs", "Ythan Estuary", "Nepean River")

# Define the colour palette for each group
site_colours <- c("#000000","#000000", "#E69F00", "#56B4E9", "#009E73",NA, "#000000", "#E69F00", "#56B4E9", "#009E73")
#and link to site names
site2colours <- cbind(site_order,site_colours)
#combine with df
df <- merge(df, site2colours,by.x="site",by.y="site_order",all.x=TRUE)
df <- unique(df)

#add null site for break
nl <- df[df$site=="Gogo Reef",]
nl$site<-""
nl$habitat <- ""
nl$group<-""
nl$value<-0
nl$site_colours<-"white"
df <- rbind(df,nl)

# Select only four specific sites to display in the legend
legend_sites <- c("Gogo Reef", "Caribbean Reefs", "Ythan Estuary", "Nepean River")

# Corresponding labels for the legend
legend_labels <- c("Devonian", "reef", "estuary","freshwater")

# Create a named vector for scale_fill_manual using all site colors
all_site_colours <- setNames(site_colours, site_order)

# Fix metric name
df$metric[df$metric == "specialization"] <- "specialisation"

# Define a scaling factor for the text size
text_scaling_factor <- 1.5

# Make light colours to use in geom_rect()
original_pink <- rgb(255/255, 182/255, 193/255)

# Mix the original pink with white to make it lighter
lighter_pink <- colorRampPalette(colors = c(original_pink, "white"))(2)[1]

# Function for plotting
plot_barchart <- function(x, metric) {
  df1 <- x[x$metric == metric, ]
  
  df1$site <- factor(df1$site, levels = site_order)
  
  # Find the maximum y value and set the y-axis limits accordingly
  max_y <- max(df1$value, na.rm = TRUE)
  y_limit <- max_y * 1.1  # Extend the limit a bit
  
  # Plotting
  ggplot(df1, aes(x = site, y = value, fill = site)) +
    # Adding transparent red rectangle
    geom_rect(aes(xmin = -Inf, xmax = "", ymin = 0, ymax = Inf), fill = rgb((255 + 255)/510, (182 + 255)/510, (193 + 255)/510)) +
    # Adding transparent yellow rectangle
    geom_rect(aes(xmin = "", xmax = Inf, ymin = 0, ymax = Inf), fill = rgb(1, 1, 0.5)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = all_site_colours, 
                      name = "Habitat/Group",  # Setting the legend title
                      breaks = legend_sites,  # Sites to include in the legend
                      labels = legend_labels) +  # Labels for the legend keys
    labs(x = "site", y = metric) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = text_scaling_factor * 16),  
          axis.text.y = element_text(size = text_scaling_factor * 16),  
          axis.title = element_text(size = text_scaling_factor * 18),  
          strip.text = element_text(size = text_scaling_factor * 18),  
          legend.text = element_text(size = text_scaling_factor * 18),  
          legend.title = element_text(size = text_scaling_factor * 20),  
          plot.title = element_text(size = text_scaling_factor * 20),   
          legend.position = c(0.825, 0.9),  
          legend.justification = c(0.5, 1),  
          legend.box.background = element_rect(color = "black", size = 0.5),  
          legend.key.size = unit(1, "cm"),
          legend.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(color = NA, size = NA),  # Add a black border around the panel
          plot.background = element_rect(color = NA, fill = "white")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, y_limit))
}

ric <- plot_barchart(df,"functional richness") + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
nn <- plot_barchart(df,"nearest neighbour") + theme(legend.position = "none",  axis.title.x = element_blank(),axis.text.x = element_blank())
spe <- plot_barchart(df,"specialisation") + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank())
div <- plot_barchart(df,"divergence") + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank())
eve <- plot_barchart(df,"evenness") + theme(legend.position = "none", axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))


# Combine the ggplots and save as pdf (ric, ric2, nn, nn2)
# Define the heights for each row
row_heights <- c(0.25,2.4,3.1)  


## Do the same for the other metrics to go in supplementary material
# Combine the ggplots and save as pdf (spe, spe2, eve, eve2, div, div2)
# Define the heights for each row
row_heights <- c(0.25,2.4,2.4,3.2)  # Adjust the heights
#add panel labels
p1 <- grid.text("c.", y = 59.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
#trop <- grid.text("tropical", y = 50, x = 0.35,gp = gpar(fontsize = 30, fontface = "bold", col = "red") )
#temp <- grid.text("temperate/\nsubtropical", y = 55, x = 0.85,gp = gpar(fontsize = 30, fontface = "bold", col = "grey22") )
speA <- grid.arrange(spe, p1, ncol = 1, heights = c(3, 0.05,0.05,0.05))
p2 <- grid.text("e.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
eveA <- grid.arrange(eve, p2, ncol = 1, heights = c(1, 0.05))  
p3 <- grid.text("d.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
divA <- grid.arrange(div, p3, ncol = 1, heights = c(1, 0.05))  

p4 <- grid.text("a.", y = 59.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
trop <- grid.text("tropical", y = 45, x = 0.475,gp = gpar(fontsize = 30, fontface = "bold", col = "red") )
temp <- grid.text("temperate/\nsubtropical", y = 30, x = 0.80,gp = gpar(fontsize = 30, fontface = "bold", col = "grey22") )
ricA <- grid.arrange(ric, p4,trop,temp, ncol = 1, heights = c(3, 0.05,0.05,0.05))  
p5 <- grid.text("b.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
nnA <- grid.arrange(nn, p5, ncol = 1, heights = c(1, 0.05))  
p6 <- grid.text("b.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
nnA <- grid.arrange(nn, p6, ncol = 1, heights = c(1, 0.05))  

bp1 <- ggplot() + 
  annotate("text", label = "observed metric", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

bp2 <- ggplot() + 
  annotate("text", label = "observed metric", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

g1 <- grid.arrange(bp1,bp2, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g2 <- grid.arrange(ricA,nnA, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g3 <- grid.arrange(speA,divA, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
empty_space <- grid::nullGrob()
g4 <- grid.arrange(eveA, eveA, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))

# Now use grid.arrange to layout the two rows
combined_plots <- grid.arrange(
  g1,
  g2,
  g3,
  g4,
  heights = row_heights,  # Use the defined row heights
  nrow = 4
)

ggsave(plot = combined_plots, 
       filename = "filepath/supplementary material/figures/figure #. Three FD metrics all species High AccuracyS.png",
       height = 32, width =20,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")
