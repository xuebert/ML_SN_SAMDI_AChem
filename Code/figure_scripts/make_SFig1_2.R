# this script calculates p-values for AA enrichment in the K-array and H-array
# this script also generates the p-value plots in SFig1 and SFig2

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())

if("plotrix" %in% rownames(installed.packages()) == FALSE) {install.packages("plotrix")}
library(plotrix)

source("support_functions/color_to_hex.R")
source("support_functions/data_to_p_values.R")

load("../Data/K_array_data.RData")
load("../Data/H_array_data.RData")
SN_storage = list(colMeans(SN_K), colMeans(SN_H))

AA_names = c("A", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") # in alphabetical order 

########################################################################################

# stat parameters
location_size = matrix(c(50, 25, 40, 25), ncol = 2)

# significance threshold
p_value = 0.0001
threshold = -log10(p_value)

# storage
storage = vector("list", 2)
AA_table_p = storage
AA_table = storage
AA_tot = storage

# calculate p values
for (n_array in 1:2) {
  AA_storage = t(sapply(1:361, function(i) strsplit(names(SN_storage[[n_array]]), "")[[i]]))
  if (n_array == 2) {AA_storage = AA_storage[,1:2]}
  
  return_list = data_to_p_values(data_value = SN_storage[[n_array]], AA_storage = AA_storage, low_threshold = location_size[[1, n_array]], high_threshold = location_size[[2, n_array]])
  
  AA_table_p[[n_array]] = return_list[[1]]
  AA_table[[n_array]] = return_list[[2]]
  AA_tot[[n_array]] = return_list[[3]]
}

########################################################################################
# plotting

# plot parameters
xlim = c(0, 8)
colors = c("blue", "red")
axis_scale = seq(0, 15, by = 5)

pdf(file = "figure_dump/SFig2_AA_enrichment", height = 10, width = 18)

# margins
par(mfrow = c(1,2))
# other plot parameters
par(cex.main = 1, mai = c(1.5, 1.2, 1, 0.5), ps = 24)
par(las = 2)

array_names = c("K", "H")
for (n_array in 1:2) {
  
  ### low quality ### 
  xlab = "-log10(p value)"
  ylab = "Amino acids"
  main = paste(array_names[[n_array]], "low signal to noise ratio")
  
  # generate barplots
  barplot(as.table(-log10(AA_table_p[[n_array]][2:1,19:1])), horiz = T, beside = T, axes = F, xlim = xlim, xlab = xlab, ylab = ylab, main = main, xaxt = "n", cex.names = 0.75, col = colors, mgp = c(4, 1, 0))
  par(las = 1)
  axis(side = 1, at = 0:8, labels = 0:8)
  lines(x = c(threshold, threshold), y = c(0, 58), col = "black", lty = 2)
  
  ### high quality ###
  ylab = ""
  main = paste(array_names[[n_array]], "high signal to noise ratio")
  
  par(cex.main = 1, mai = c(1.5, 0.5, 1, 1.2), ps = 24)
  barplot(as.table(-log10(AA_table_p[[n_array]][4:3,19:1])), horiz = T, beside = T, axes = F, xlim = xlim, xlab = xlab, ylab = ylab, main = main, xaxt = "n", cex.names = 0.75, col = colors, mgp = c(4, 1, 0))
  
  axis(side = 1, at = 0:8, labels = 0:8)
  lines(x = c(threshold, threshold), y = c(0, 58), col = "black", lty = 2)
}
dev.off()


