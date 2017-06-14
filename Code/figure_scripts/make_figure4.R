# this script makes the barplots found in Fig4 and SFig4

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())

load("../Data/K_array_data.RData")
load("../Data/H_array_data.RData")
SN_storage = list(colMeans(SN_K), colMeans(SN_H))
AA_storage = list(t(sapply(1:361, function(i) strsplit(names(SN_storage[[1]]), "")[[i]])), t(sapply(1:361, function(i) strsplit(names(SN_storage[[2]]), "")[[i]])))
AA_storage[[2]] = AA_storage[[2]][,1:2]

# rename peptides in 2nd array (get rid of trailing H)
names(SN_storage[[2]]) = sapply(1:361, function(i) paste(AA_storage[[2]][i,], collapse = ""))

# rearrange SN values into grid
AA_names = c("A", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (n_array in 1:2) {
  # sort into alphabetical order
  SN_storage[[n_array]] = SN_storage[[n_array]][order(names(SN_storage[[n_array]]))]
  
  SN_storage[[n_array]] = matrix(SN_storage[[n_array]], ncol = 19, byrow = T)
  SN_storage[[n_array]] = SN_storage[[n_array]] / max(SN_storage[[n_array]])
  
  # rename columns and rows into alphabetical order
  rownames(SN_storage[[n_array]]) = AA_names
  colnames(SN_storage[[n_array]]) = AA_names
  
}

###########################################
# bubble chart

xlab = "Signal to noise ratio"
ylab = "Amino acids"
main = ""
colors = "black"
xlim = c(0, 1.05)
array_names = c("Kac", "K")

AA_order = order(colMeans(SN_storage[[1]]) + rowMeans(SN_storage[[1]]))

pdf("figure_dump/Fig4_SFig4_barplots.pdf")
for (n_array in 1:2) {
  
  plot_storage = (colMeans(SN_storage[[n_array]]) + rowMeans(SN_storage[[n_array]]))[AA_order]
  
  main = paste(array_names[[n_array]], "XZ-positions")
  
  # barplot for both positions together
  barplot(as.table(rev(plot_storage)), horiz = T, beside = T, axes = F, xlim = xlim, xlab = xlab, ylab = ylab, main = main, xaxt = "n", cex.names = 0.75, col = colors, mgp = c(3, 1, 0), space = 5)
  box(col = "black")
  
  # also make it for each x and z positions.
  # first x position
  plot_storage = colMeans(SN_storage[[n_array]])[AA_order]
  plot_storage = plot_storage / max(plot_storage)
  
  main = paste(array_names[[n_array]], "X-position")
  barplot(as.table(rev(plot_storage)), horiz = T, beside = T, axes = F, xlim = xlim, xlab = xlab, ylab = ylab, main = main, xaxt = "n", cex.names = 0.75, col = colors, mgp = c(3, 1, 0), space = 5)
  box(col = "black")
  
  # z position
  plot_storage = rowMeans(SN_storage[[n_array]])[AA_order]
  plot_storage = plot_storage / max(plot_storage)
  
  main = paste(array_names[[n_array]], "Z-position")
  barplot(as.table(rev(plot_storage)), horiz = T, beside = T, axes = F, xlim = xlim, xlab = xlab, ylab = ylab, main = main, xaxt = "n", cex.names = 0.75, col = colors, mgp = c(3, 1, 0), space = 5)
  box(col = "black")
}

dev.off()
