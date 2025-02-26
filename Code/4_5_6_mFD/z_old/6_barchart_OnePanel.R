#barchart of mFD results
library(mFD)
library(dplyr)
library(tidyr)
library(gawdis)
library(ggpubr)
library(grid)
library(gridExtra)
library(patchwork)
library(colorBlindness)

#Adjust file path at"###"

setwd("###")
Plot1 <- readRDS("###/data/data_for_plots/Plot1.RDS")
Plot3 <- readRDS("###/data/data_for_plots/Plot3.RDS")

#cut to three metrics and make Devonian black so they stand out, tropical with a black outline
d1 <- Plot1$data
d3 <- Plot3$data
d1 <- d1[d1$metric!="originality",]
d3 <- d3[d3$metric!="originality",]

#assign data 
df <- Plot3$data
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

# Define the order of sites within each group
site_order <- c("Gogo Reef","Miguasha Estuary", "Caribbean Reefs", "Santa Cruz Estuary", "Braço Morto\nAcima and Abaixo","","Chile Reefs", "Ythan Estuary", "Nepean River")

# Define the colour palette for each group
site_colours <- c("#000000","#000000", "#E69F00", "#56B4E9", "#009E73",NA, "#E69F00", "#56B4E9", "#009E73")
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

#fix metric name
df$metric[df$metric == "specialization"] <- "specialisation"

# Define a scaling factor for the text size
text_scaling_factor <- 1.5

#make light colours to use in geom_rect()
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
nn <- plot_barchart(df,"nearest neighbour") + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank())
spec <- plot_barchart(df,"specialisation") + theme(legend.position = "none", axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))

#now do it for metrics that were calculated controlling for species diversity
df2 <- Plot1$data
df2$metric <- as.character(df2$metric)

df2$metric <- ifelse(df2$metric=="richness", "functional richness", 
                    ifelse(df2$metric=="nearest\nneighbour", "nearest neighbour",df2$metric))

df2$site <- ifelse(df2$site == "BracoMorto", "Braço Morto\nAcima and Abaixo",
                  ifelse(df2$site == "Caribbean", "Caribbean Reefs",
                         ifelse(df2$site == "Chile_reef", "Chile Reefs", 
                                ifelse(df2$site == "Nepean", "Nepean River",
                                       ifelse(df2$site == "Santa_Cruz_Channel", "Santa Cruz Estuary",
                                              ifelse(df2$site == "Ythan","Ythan Estuary", df2$site))))))

df2$site[df2$site=="Gogo"] <- "Gogo Reef"
df2$site[df2$site=="Miguasha"] <- "Miguasha Estuary"

df2$value[df2$metric=="functional richness"] <- df2$value[df2$metric=="functional richness"]*100000

#combine colours with df2
df2 <- merge(df2, site2colours,by.x="site",by.y="site_order",all.x=TRUE)
df2 <- unique(df2)

#add null site for break
df2 <- rbind(df2,nl)

#fix metric name
df2$metric[df2$metric == "specialization"] <- "specialisation"

ric2 <- plot_barchart(df2,"functional richness") + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_blank())
nn2 <- plot_barchart(df2,"nearest neighbour") + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank(),axis.title.y=element_blank())
spec2 <- plot_barchart(df2,"specialisation") + theme(legend.position = "none",axis.title.y=element_blank(),axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))

#combine the ggplots and save as pdf (ric, ric2, nn, nn2, sp, sp2)
# Define the heights for each row
row_heights <- c(0.25,2.4,2.4,3.2)  # Adjust the heights
#add panel labels
p1 <- grid.text("a.", y = 59.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
trop <- grid.text("tropical", y = 45, x = 0.5,gp = gpar(fontsize = 30, fontface = "bold", col = "red") )
temp <- grid.text("temperate/\nsubtropical", y = 30, x = 0.80,gp = gpar(fontsize = 30, fontface = "bold", col = "grey22") )
ricA <- grid.arrange(ric, p1,trop,temp, ncol = 1, heights = c(3, 0.05,0.05,0.05))
p2 <- grid.text("c.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
nnA <- grid.arrange(nn, p2, ncol = 1, heights = c(1, 0.05))  
p3 <- grid.text("e.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
specA <- grid.arrange(spec, p3, ncol = 1, heights = c(1, 0.05))  

p4 <- grid.text("b.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
ric2A <- grid.arrange(ric2, p4, ncol = 1, heights = c(1, 0.05))  
p5 <- grid.text("d.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
nn2A <- grid.arrange(nn2, p5, ncol = 1, heights = c(1, 0.05))  
p6 <- grid.text("f.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
spec2A <- grid.arrange(spec2, p6, ncol = 1, heights = c(1, 0.05))  

bp1 <- ggplot() + 
  annotate("text", label = "all species", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

bp2 <- ggplot() + 
  annotate("text", label = "diversity-controlled", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

g1 <- grid.arrange(bp1,bp2, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g2 <- grid.arrange(ricA,ric2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g3 <- grid.arrange(nnA, nn2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g4 <- grid.arrange(specA, spec2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))

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
       filename = "###/figure 2. Three FD metrics all species and control for diversity.png",
       height = 32, width =20,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white"
)

