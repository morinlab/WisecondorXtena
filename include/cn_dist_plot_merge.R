#!/usr/bin/env Rscript

library(mclust)
library(tidyverse)
#library(glue)
library(ggplot2)
library(cowplot)
library(png)
library(grid)
library(patchwork)

args <-  commandArgs(trailingOnly = TRUE)

# Test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Input bins.bed file and PNG must be supplied (input file).n", call.=FALSE)
}

# Import bins.bed and genome_wide.png files
in_bed <- read.delim(args[1])
CN_plot <- readPNG(args[2])

# Get rid of all rows that have blacklisted values (called NaN)
in_bed_nonull <- na.omit(in_bed, 
                         cols_only(ratio))

# Fit mixture model to the log2 ratio portion of the cleaned up input bed file
model_fit = Mclust(in_bed_nonull$ratio, 
                           G=1, 
                           warn = TRUE)

# Extract the mean of the classification/mixture model and make it the variable called "modeldf"
#model_fit$parameters$mean
model_mean=data.frame(mean=model_fit$parameters$mean)

# Plot the distribution curve and rotate horizontally
dist_plot <- 
  ggplot(in_bed_nonull, aes(x=ratio)) + 
  geom_density() + 
  geom_vline(aes(xintercept = 0, 
                 colour = "Neutral"), 
             linetype="dotted", 
             size = 0.35, 
             show.legend = TRUE, 
             key_glyph="path") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = c(0.8,0.8), 
        legend.title = element_blank()) +
  geom_vline(data=model_mean, 
             aes(xintercept=mean, 
                 colour = "Offset"), 
             linetype="longdash", 
             size = 0.5,
             show.legend = TRUE, 
             key_glyph = "path") +
  xlim(-2.2,2.7) +
  #ylim(0,2) +
  scale_color_manual(values = c(Offset = "black", Neutral = "blue")) +
  coord_flip()


# Give copy number plot a labeled y-axis that matches up to the same y-axis as the density plot, uses the grid package
CN_axis_align <- ggplot(in_bed_nonull, aes(x=ratio)) +
  xlim(0, 10) +
  ylim(-2.2, 2.7) +
  theme(axis.line = element_line(),
        axis.line.x = element_line(color = "white"),
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),
        axis.ticks.x = element_line(color = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(y = "log2 ratio") + 
  annotation_custom(rasterGrob(CN_plot, width = 1, height = 1),
                    xmin = -1.2, xmax = 10.8,
                    ymin = -6.07, ymax = 2.87)

# Remove y-axis labels from distribution plot
dist_plot_no_y <- 
  dist_plot +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"))

# Align both plots together using the patchwork package
pdf(args[3], width = 8, height = 4)

plot(CN_axis_align + 
  dist_plot_no_y + 
  plot_layout(widths = c(4.5,1)))

dev.off()
