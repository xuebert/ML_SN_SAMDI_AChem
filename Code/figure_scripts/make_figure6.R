# this script combines and makes the ribbon plots for Q2 values when peptides are sampled at different sizes (figure6)

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())

source("support_functions/color_to_hex.R")

# grab property names
source("../Data/amino_acid_physical_properties/AA_physical_properties.R")
property_names = names(AA_physical_properties(AA_list = "A", full_name = T))

# load results
load("sampling_peptides_results.RData")

color1 = color_to_hex("red", transparency = 0.7)
color2 = color_to_hex("blue", transparency = 0.7)
color3 = color_to_hex("grey80", transparency = 0.3)

########################################
# 
pdf("figure_dump/Fig6_peptides_sampled.pdf")
xlim = c(0, 361)
ylim = c(-0.2, 1)
sample_range = c(5:350)
LB = 3 # 10% confidence interval
UB = 18 # 90% confidence interval (total of 80% confidence interval around the median value)

make_plot <- function(n_array = 1, n_random = 1, n_color) {
  
  plot_storage = data.frame(num_samples = sample_range, median = Q2_sampled[[n_random]][[n_array]][11, sample_range], min =  Q2_sampled[[n_random]][[n_array]][LB, sample_range], max = Q2_sampled[[n_random]][[n_array]][UB, sample_range])
  
  polygon(x = c(rev(sample_range), sample_range), y = c(rev(plot_storage$min), plot_storage$max), col = n_color, border = NA)
  par(new = T)
  plot(x = sample_range, plot_storage$median, type = "l", xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", col = substr(n_color, 1, 7))
  
}

# actual K-array data
plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = "Number of sampled peptides", ylab = "Q2", main = "Predicting with limited peptide samples", xaxt = "n")
par(new = T)
make_plot(n_color = color1)
axis(side = 1, at = c(0, 100, 200, 300, 361), labels = c(0, 100, 200, 300, 361))

# actual H-array data
par(new = T)
make_plot(n_array = 2, n_color = color2)

# random data
par(new = T)
make_plot(n_array = 1, n_random = 2, n_color = color3)

lines(x = c(-20, 380), y = c(0, 0), lty = 2)

# calculate shortest distances from c(1,1)
distances_K = ((1 - Q2_sampled$actual$K[11,])^2 + ((0 - 1:361) / 361) ^ 2) ^ 0.5
optimal_point_K = order(distances_K)[[1]]
optimal_point_Q2_K = Q2_sampled[[1]][[1]][[11, optimal_point_K]]
arrows(x0 = optimal_point_K, x1 = optimal_point_K, y0 = optimal_point_Q2_K + 0.1, y1 = optimal_point_Q2_K, length = 0.2)

distances_H = ((1 - Q2_sampled$actual$H[11,])^2 + ((0 - 1:361) / 361) ^ 2) ^ 0.5
optimal_point_H = order(distances_H)[[1]]
optimal_point_Q2_H = Q2_sampled[[1]][[2]][[11, optimal_point_H]]
arrows(x0 = optimal_point_H, x1 = optimal_point_H, y0 = optimal_point_Q2_H - 0.1, y1 = optimal_point_Q2_H - 0.02, length = 0.2)

legend("topleft", legend = c("K", "H", "Random"), fill = c(color1, color2, color3))

dev.off()

