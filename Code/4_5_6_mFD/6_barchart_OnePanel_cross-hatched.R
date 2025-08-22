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
library(ggpattern)

#Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)

Plot3 <- readRDS("filepath/output/Plot3.RDS")
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
site_colours <- c("#E69F00","#56B4E9", "#E69F00", "#56B4E9", "#009E73",NA, "#009E73", "#E69F00", "#56B4E9", "#009E73")
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
legend_sites <- c("Caribbean Reefs", "Ythan Estuary", "Nepean River")

# Corresponding labels for the legend
legend_labels <- c("reef", "estuary","freshwater")

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

# assign which get hatching
hatch_sites <- c("Gogo Reef", "Miguasha Estuary", "Canowindra Billabong")

# Function for plotting
plot_barchart <- function(x, metric) {
  df1 <- x[x$metric == metric, ]
  df1$site <- factor(df1$site, levels = site_order)
  
  # pattern flag: "stripe" for selected bars, "none" otherwise
  df1$pattern <- ifelse(df1$site %in% hatch_sites, "stripe", "none")
  df1$pattern <- factor(df1$pattern, levels = c("none","stripe"))
  
  max_y <- max(df1$value, na.rm = TRUE)
  y_limit <- max_y * 1.1
  
  ggplot(df1, aes(x = site, y = value, fill = site)) +
    # background rectangles (unchanged)
    geom_rect(aes(xmin = -Inf, xmax = "", ymin = 0, ymax = Inf),
              fill = rgb((255 + 255)/510, (182 + 255)/510, (193 + 255)/510)) +
    geom_rect(aes(xmin = "", xmax = Inf, ymin = 0, ymax = Inf),
              fill = rgb(1, 1, 0.5)) +
    geom_hline(yintercept = 0, color = "black") +
    
    # patterned bars instead of geom_bar
    ggpattern::geom_col_pattern(
      aes(pattern = pattern),
      position = "dodge",
      # pattern appearance (tweak to taste)
      pattern_angle = 45,          # diagonal hatching
      pattern_spacing = 0.03,      # distance between stripes (in npc units)
      pattern_density = 0.5,       # stripe thickness as fraction of spacing
      pattern_fill = "black",      # colour of the hatch lines
      pattern_alpha = 1,
      pattern_key_scale_factor = 0.6
    ) +
    
    scale_fill_manual(values = all_site_colours,
                      name = "habitat",
                      breaks = legend_sites,
                      labels = legend_labels) +
    
    # separate legend for hatch pattern (Devonian vs Modern)
    scale_pattern_manual(
      name = "period",
      values = c("none" = "none", "stripe" = "stripe"),
      labels = c("modern", "Devonian")
    ) +
    
    labs(x = "site", y = metric) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = text_scaling_factor * 16),
      axis.text.y = element_text(size = text_scaling_factor * 16),
      axis.title  = element_text(size = text_scaling_factor * 18),
      strip.text  = element_text(size = text_scaling_factor * 18),
      legend.text = element_text(size = text_scaling_factor * 18),
      legend.title= element_text(size = text_scaling_factor * 20),
      plot.title  = element_text(size = text_scaling_factor * 20),
      
      # --- Option A: one unified box around both legends ---
      legend.box = "vertical",                    # stack the two guides
      legend.box.just = "left",
      legend.background = element_blank(),        # remove per-guide boxes
      legend.box.background = element_rect(       # single outer box
        colour = "black", fill = "white", linewidth = 0.5
      ),
      legend.box.margin = margin(6, 6, 6, 6),
      
      # consistent key sizing
      legend.key.size  = unit(1, "cm"),
      legend.key.width = unit(1, "cm"),
      
      # place the unified legend box
      legend.position = c(0.825, 0.975),
      legend.justification = c(0.5, 1),
      
      panel.background = element_rect(color = NA, size = NA),
      plot.background  = element_rect(color = NA, fill = "white")
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, y_limit)) +
    
    guides(
      # Habitat legend: solid colours (no hatch)
      fill = guide_legend(
        order = 1,
        ncol = 1,
        keywidth = unit(1, "cm"),
        title.position = "top",
        override.aes = list(pattern = "none")
      ),
      # Pattern legend: hatch swatches
      pattern = guide_legend(
        order = 2,
        ncol = 1,
        keywidth = unit(1, "cm"),
        title.position = "top",
        override.aes = list(
          fill = "grey85",
          pattern_fill = "black",
          pattern_background = NA
        )
      )
    )
}


ric <- plot_barchart(df,"functional richness") + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
nn <- plot_barchart(df,"nearest neighbour") + theme(legend.position = "none", axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))
spe <- plot_barchart(df,"specialisation") + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_blank())
div <- plot_barchart(df,"divergence") + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank())
eve <- plot_barchart(df,"evenness") + theme(legend.position = "none", axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))

#SES data
Plot1 <- readRDS("filepath/output/Plot1_SES.RDS")
d1 <- unique(Plot1$data)

# SES by climate zone
src <- d1[d1$metric=="richness",]
mean(src$value[src$site%in%c("Canowindra","Chile_reef","Nepean","Ythan")])
mean(src$value[!src$site%in%c("Canowindra","Chile_reef","Nepean","Ythan")])


#now do it for metrics that were calculated controlling for species diversity
df2 <- d1 
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
df2$site[df2$site=="Canowindra"] <- "Canowindra Billabong"


#combine colours with df2
df2 <- merge(df2, site2colours,by.x="site",by.y="site_order",all.x=TRUE)
df2 <- unique(df2)

#add null site for break
df2 <- rbind(df2,nl)

#fix metric name
df2$metric[df2$metric == "specialization"] <- "specialisation"
df2$metric[df2$metric == "eveness"] <- "evenness"

# Function for plotting when values could be positive, negative, or mixed
plot_barchart2 <- function(x, metric) {
  # Subset data for the selected metric
  df1 <- x[x$metric == metric, ]
  
  # Ensure site is a factor with your specified order
  df1$site <- factor(df1$site, levels = site_order)
  
  # Create a clean pattern flag; avoid NAs (e.g., spacer site == "")
  df1$pattern <- ifelse(df1$site %in% hatch_sites & df1$site != "" & !is.na(df1$value),
                        "stripe", "none")
  df1$pattern <- factor(df1$pattern, levels = c("none", "stripe"))
  df1 <- droplevels(df1)
  
  # y-limits with buffer
  min_y <- min(df1$value, na.rm = TRUE)
  max_y <- max(df1$value, na.rm = TRUE)
  buffer <- (max_y - min_y) * 0.1
  y_limit <- c(min_y - buffer, max_y + buffer)
  
  ggplot(df1, aes(x = site, y = value, fill = site)) +
    # Background halves (adjust xmin/xmax to your split index)
    annotate("rect", xmin = 6, xmax = Inf, ymin = min_y - buffer, ymax = max_y + buffer,
             fill = rgb(1, 1, 0.5)) +
    annotate("rect", xmin = -Inf, xmax = 6, ymin = min_y - buffer, ymax = max_y + buffer,
             fill = rgb((255 + 255)/510, (182 + 255)/510, (193 + 255)/510)) +
    geom_hline(yintercept = 0, color = "black") +
    
    # PATTERNED bars (replaces geom_bar)
    ggpattern::geom_col_pattern(
      aes(pattern = pattern),
      position = "dodge",
      # pattern appearance (tweak to taste)
      pattern_angle = 45,          # diagonal hatching
      pattern_spacing = 0.03,      # distance between stripes (in npc units)
      pattern_density = 0.5,       # stripe thickness as fraction of spacing
      pattern_fill = "black",      # colour of the hatch lines
      pattern_alpha = 1,
      pattern_key_scale_factor = 0.6
    ) +
    
    # Colours for Habitat/Group
    scale_fill_manual(values = all_site_colours,
                      name = "habitat",
                      breaks = legend_sites,
                      labels = legend_labels) +
    
    # Pattern legend (time period)
    scale_pattern_manual(
      name = "time period",
      values = c("none" = "none", "stripe" = "stripe"),
      labels = c("Modern", "Devonian")
    ) +
    
    labs(x = "Site", y = metric) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = y_limit) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = text_scaling_factor * 16),
      axis.text.y = element_text(size = text_scaling_factor * 16),
      axis.title   = element_text(size = text_scaling_factor * 18),
      legend.text  = element_text(size = text_scaling_factor * 18),
      legend.title = element_text(size = text_scaling_factor * 20),
      legend.position = c(0.825, 0.9),
      legend.justification = c(0.5, 1),
      legend.key.size = unit(1, "cm"),
      legend.background = element_rect(fill = "white", colour = "white"),
      legend.box.background = element_rect(colour = "black", size = 0.5) # border around legend
    ) +
    # Make Habitat/Group legend show solid colours (no hatch);
    # keep a separate pattern legend for time period
    guides(
      fill = guide_legend(override.aes = list(pattern = "none"), order = 1),
      pattern = guide_legend(
        override.aes = list(fill = "grey85", pattern_colour = "black"),
        order = 2
      )
    )
}


# Make the plot objects
ric2 <- plot_barchart2(df2,"functional richness") + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank())
nn2 <- plot_barchart2(df2,"nearest neighbour") + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))
spe2 <- plot_barchart2(df2,"specialisation") + theme(legend.position = c(0.4,0.9), axis.title.x = element_blank(),axis.text.x = element_blank())
div2 <- plot_barchart2(df2,"divergence") + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank())
eve2 <- plot_barchart2(df2,"evenness") + theme(legend.position = "none", axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1))

# Combine the ggplots and save as pdf (ric, ric2, nn, nn2)
# Define the heights for each row
row_heights <- c(0.25,2.4,3.1)  
# Add panel labels
p1 <- grid.text("a.", y = 59.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
trop <- grid.text("tropical", y = 45, x = 0.475,gp = gpar(fontsize = 30, fontface = "bold", col = "red") )
temp <- grid.text("temperate/\nsubtropical", y = 30, x = 0.80,gp = gpar(fontsize = 30, fontface = "bold", col = "grey22") )
ricA <- grid.arrange(ric, p1,trop,temp, ncol = 1, heights = c(3, 0.05,0.05,0.05))
p2 <- grid.text("c.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
nnA <- grid.arrange(nn, p2, ncol = 1, heights = c(1, 0.05))  

p4 <- grid.text("b.", y = 20.5, x = 0.1,gp = gpar(fontsize = 30, fontface = "bold") )
ric2A <- grid.arrange(ric2, p4, ncol = 1, heights = c(1, 0.05))  
p5 <- grid.text("d.", y = 20.5, x = 0.1,gp = gpar(fontsize = 30, fontface = "bold") )
nn2A <- grid.arrange(nn2, p5, ncol = 1, heights = c(1, 0.05))  

bp1 <- ggplot() + 
  annotate("text", label = "observed metric", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

bp2 <- ggplot() + 
  annotate("text", label = "standardized effect size", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

g1 <- grid.arrange(bp1,bp2, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g2 <- grid.arrange(ricA,ric2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g3 <- grid.arrange(nnA, nn2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))

# Now use grid.arrange to layout the two rows
combined_plots <- grid.arrange(
  g1,
  g2,
  g3,
  heights = row_heights,  # Use the defined row heights
  nrow = 3
)

ggsave(plot = combined_plots, 
       filename = "filepath/figures/for manuscript/figure 1. Two FD metrics all species and control for diversity_wSES.png",
       height = 22, width =20,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white"
)

## Do the same for the other metrics to go in supplementary material
# Combine the ggplots and save as pdf (spe, spe2, eve, eve2, div, div2)
# Define the heights for each row
row_heights <- c(0.25,2.4,2.4,3.2)  # Adjust the heights
#add panel labels
p1 <- grid.text("a.", y = 59.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
trop <- grid.text("tropical", y = 50, x = 0.35,gp = gpar(fontsize = 30, fontface = "bold", col = "red") )
temp <- grid.text("temperate/\nsubtropical", y = 55, x = 0.85,gp = gpar(fontsize = 30, fontface = "bold", col = "grey22") )
speA <- grid.arrange(spe, p1,trop,temp, ncol = 1, heights = c(3, 0.05,0.05,0.05))
p2 <- grid.text("e.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
eveA <- grid.arrange(eve, p2, ncol = 1, heights = c(1, 0.05))  
p3 <- grid.text("c.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
divA <- grid.arrange(div, p3, ncol = 1, heights = c(1, 0.05))  

p4 <- grid.text("b.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
spe2A <- grid.arrange(spe2, p4, ncol = 1, heights = c(1, 0.05))  
p5 <- grid.text("f.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
eve2A <- grid.arrange(eve2, p5, ncol = 1, heights = c(1, 0.05))  
p6 <- grid.text("d.", y = 20.5, x = 0.15,gp = gpar(fontsize = 30, fontface = "bold") )
div2A <- grid.arrange(div2, p6, ncol = 1, heights = c(1, 0.05))  

bp1 <- ggplot() + 
  annotate("text", label = "observed metric", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

bp2 <- ggplot() + 
  annotate("text", label = "standardised effect size", x = 0.5, y = 0.6, angle = 0, size = 14) +
  xlim(0, 1) + 
  ylim(0.575, 0.612) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

g1 <- grid.arrange(bp1,bp2, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g2 <- grid.arrange(speA,spe2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g3 <- grid.arrange(divA, div2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))
g4 <- grid.arrange(eveA, eve2A, ncol = 3, widths = c(3, 0.3, 2.9), layout_matrix = rbind(c(1, NA, 2)))

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
       filename = "filepath/supplementary material/figures/Sfigure 6. Three FD metrics all species and control for diversity_SES.png",
       height = 32, width =20,  dpi = 300, device = "png",limitsize = FALSE,
       bg = "white")