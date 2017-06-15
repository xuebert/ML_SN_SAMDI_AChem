# make figure 2 SN_distribution

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())

library(plotrix)

source("support_functions/color_to_hex.R")

load("../Data/K_array_data.RData")
load("../Data/H_array_data.RData")
SN_storage = list(colMeans(SN_K), colMeans(SN_H))
AA_storage = list(t(sapply(1:361, function(i) strsplit(names(SN_storage[[1]]), "")[[i]])), t(sapply(1:361, function(i) strsplit(names(SN_storage[[2]]), "")[[i]])))
AA_storage[[2]] = AA_storage[[2]][,1:2]

################### find low/high SN region size #########################

# put in sliding window slope calculation to find
window_size = 15

# storage for the different lm objects
lm_storage = vector("list", 361 - window_size)
slope_storage = vector("numeric", 361 - window_size)

for (m_index in (window_size+1):361) {
  X = (m_index-window_size):m_index
  Y = sort(SN_storage[[2]])[X]
  lm_storage[[m_index]] = lm(Y~X)
  slope_storage[[m_index]] = lm_storage[[m_index]]$coefficients[[2]]
}

plot(slope_storage) # look for when slope changes sharply

find_indices <- function(n_array = 1, n_position = 1, AA) {
  return(which(AA_storage[[n_array]][,n_position] == AA))
}

########################################################################
# # plot the distribution
pdf(file = "figure_dump/Fig2_SN_distribution.pdf")

location_size = matrix(c(50, 25, 40, 25), ncol = 2) # found from plotting slopes
ylim = list(c(0, 325), c(0,275))
axis_storage = list(list(side=2, at=c(0,100,200,300,325), labels=c(0,100,200,300,325)), list(side=2,at=c(0,100,200,275), labels=c(0,100,200,275)))

color_list = list(c("red", "red", "green", "blue"), c("red", "red", "red", "green", "blue"))
for (n_array in 1:2) {
  
  # find the three amino acid positions to serve as control.  Two are the p values closest to the threshold.  
  # X-M is closest to low SN.  X-D for random middle, Z-T for closest to high SN
  if (n_array == 1) {
    significant_AAs = list(find_indices(n_array, 1, "W"), find_indices(n_array, 1, "L"), find_indices(n_array, 2, "G"), find_indices(n_array, 1, "A"))
    significant_points = unique(unlist(significant_AAs)) 
  } else {
    significant_AAs = list(find_indices(n_array, 2, "D"), find_indices(n_array, 1, "E"), find_indices(n_array, 1, "D"), find_indices(n_array, 2, "F"), find_indices(n_array, 1, "A"))
    significant_points = unique(unlist(significant_AAs)) 
  }
  
  SN_order = order(order(SN_storage[[n_array]], decreasing = F), decreasing = F)
  
  plot(x = SN_order[-significant_points], y = SN_storage[[n_array]][-significant_points], xlab = "Sorted peptides", ylab = "Mean signal to noise ratio", main = "Signal to noise ratio distribution", xlim = c(0, 380), ylim = ylim[[n_array]], pch = "o", cex = 0.5, col = "black", yaxt = "n")
  do.call(axis, axis_storage[[n_array]])
  
  for (n_sign_AAs in 1:length(color_list[[n_array]])) {
    text(x = SN_order[significant_AAs[[n_sign_AAs]]], y = SN_storage[[n_array]][significant_AAs[[n_sign_AAs]]], labels = "o", col = color_list[[n_array]][[n_sign_AAs]], cex = 1.1)
  }
  
  gray_color = color_to_hex("grey", 0.4)
  
  # make fitted plot for low and high regions
  X = 1:location_size[[1, n_array]]
  Y = sort(SN_storage[[n_array]], decreasing = F)[X]
  lm_plot = lm(Y~X)
  x_plot_values = X[c(1, length(X))]
  y_plot_values = lm_plot$coefficients[[1]] + lm_plot$coefficients[[2]] * x_plot_values
  lines(y = y_plot_values, x = x_plot_values, col = gray_color, lwd = 20)
  
  X = 361:(361-location_size[[2, n_array]])
  Y = sort(SN_storage[[n_array]], decreasing = F)[X]
  lm_plot = lm(Y~X)
  x_plot_values = X[c(1, length(X))]
  y_plot_values = lm_plot$coefficients[[1]] + lm_plot$coefficients[[2]] * x_plot_values
  lines(y = y_plot_values, x = x_plot_values, col = gray_color, lwd = 20)
  
}
dev.off()
