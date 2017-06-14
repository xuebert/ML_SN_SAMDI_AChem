# This script creates the violin dot plots for specific amino acids that make up figure 2

setwd("C:/Users/al/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())
assign("last.warning", NULL, envir = baseenv())

if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid")}
library(grid)
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
library(gridExtra)
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
library(ggplot2)

# gather SN values
load("../Data/K_array_data.RData")
load("../Data/H_array_data.RData")
SN_storage = list(colMeans(SN_K), colMeans(SN_H))
# store amino acid information
AA_storage = list(t(sapply(1:361, function(i) strsplit(names(SN_storage[[1]]), "")[[i]])), t(sapply(1:361, function(i) strsplit(names(SN_storage[[2]]), "")[[i]])))
AA_storage[[2]] = AA_storage[[2]][,1:2]

####
make_violin <- function(n_array = 1, n_position, AA, ymax = 325, dotsize) {
  AA_indices = which(AA_storage[[n_array]][,n_position] == AA)
  plot_storage = matrix(c(SN_storage[[n_array]][AA_indices], rep(AA, 19)), ncol = 2)
  colnames(plot_storage) = c("Numbers", "Labels")
  rownames(plot_storage) = 1:nrow(plot_storage)
  plot_storage = as.data.frame(plot_storage)
  # convert back to numeric
  plot_storage$Numbers = as.numeric(levels(plot_storage$Numbers))[plot_storage$Numbers]
  
  return_p = ggplot(as.data.frame(plot_storage), aes(Labels, Numbers)) +
    ylim(0, ymax) +
    theme(panel.background=element_blank()) +
    theme(panel.grid=element_blank()) +
    theme(legend.position="none") +
    scale_x_discrete(breaks=NULL) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background = element_rect(colour = "black")) +
    annotate("text", x = 0.6, y = 200, label = paste(n_position, AA)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = dotsize, binwidth = 20, position = "dodge", fill = "white") +
    theme(axis.ticks=element_blank())
  
  return(return_p)
}

#####################################

make_plots <- function(n_array, position_list, AA_list, filename, ymax) {
  # plot the whole distribution first
  plot_storage = matrix(SN_storage[[n_array]], ncol = 1)
  plot_storage = cbind(plot_storage, rep("SN_all", nrow(plot_storage)))
  
  colnames(plot_storage) = c("SN_all", "Labels")
  rownames(plot_storage) = 1:nrow(plot_storage)
  plot_storage = as.data.frame(plot_storage)
  # convert back to numeric
  plot_storage$SN_all = as.numeric(levels(plot_storage$SN_all))[plot_storage$SN_all]
  
  dotsize = 0.15
  p = vector("list", 5)
  p[[1]] <- ggplot(as.data.frame(plot_storage), aes(x = Labels, y = SN_all)) +
    #   geom_violin(aes(fill = factor(Labels)), trim = T) +
    ylim(0, ymax) +
    theme(panel.background=element_blank()) +
    theme(panel.grid=element_blank()) +
    theme(legend.position="none") +
    scale_x_discrete(breaks=NULL) +
    theme(axis.title.x=element_blank()) +
    theme(panel.background = element_rect(colour = "black")) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = dotsize, binwidth = 20, position = "dodge", fill = "white")
  
  # significant AAs
  for (n_plot in 2:(length(position_list) + 1)) {
    p[[n_plot]] <- make_violin(n_array = n_array, n_position = position_list[[n_plot-1]], AA = AA_list[[n_plot-1]], ymax = ymax, dotsize)
  }
  
  ggsave(filename, width = 8, height = 11.5, do.call(marrangeGrob, list(grobs = p, nrow = 1, ncol = 1)), useDingbats = FALSE)
  
}

# plot parameters
position_list = list(c(1,1,2,1), c(2,1,1,2,1)) # x or z position
AA_list = list(c("W", "L", "G", "A"), c("D", "E", "D", "F", "A")) # amino acids to plot
ymax = c(325, 275)

n_array = 1
filename = "figure_dump/Fig2_SN_violin_K.pdf"
make_plots(n_array, position_list[[n_array]], AA_list[[n_array]], filename, ymax[[n_array]])
n_array = 2
filename = "figure_dump/Fig2_SN_violin_H.pdf"
make_plots(n_array, position_list[[n_array]], AA_list[[n_array]], filename, ymax[[n_array]])
