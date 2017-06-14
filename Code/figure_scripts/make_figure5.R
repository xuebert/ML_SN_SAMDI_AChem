# this script processes the results from the feature selection
# this script creates the Q2 plots for both K- and H-array

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

rm(list=ls())

source("support_functions/color_to_hex.R")

# grab property names
source("../Data/amino_acid_physical_properties/AA_physical_properties.R")
property_names = names(AA_physical_properties(AA_list = "A", full_name = T))
property_classes = c(rep("Hydrophobic", 8), rep("Steric", 16), rep("Electronic", 15))

# load results
load("predict_SN_results.RData")
# Q2_int holds internal k-folds results within arrays.  List rows are K-array, H-array, and both arrays.
# list columns are (1) all features at once and (2) individual features
# first list column has three values: x-position, z-position, and both positions
# second list column has 1:39 for x-position, 40:78 for x-position, and 79:117 for both positions
# Q2_ext holds external validation results between arrays

################################################################################
# print out the tables of properties and Q2s (STables 1 and 2)
# K
table_storage = cbind(property_classes, property_names, Q2_int[[1,2]][1:39], Q2_int[[1,2]][40:78], Q2_int[[1,2]][79:117])
# sort by Q2 of both positions
table_storage = table_storage[order(table_storage[,5], decreasing = T), ]
table_storage = rbind(c("", "", "X-position", "Z-position", "Both positions"), table_storage)
write.table(table_storage, file = "STable1_K_Q2.csv", sep = ",", quote = F, row.names = F, col.names = F)

# H
table_storage = cbind(property_classes, property_names, Q2_int[[2,2]][1:39], Q2_int[[2,2]][40:78], Q2_int[[2,2]][79:117])
# sort by Q2 of both positions
table_storage = table_storage[order(table_storage[,5], decreasing = T), ]
table_storage = rbind(c("", "", "X-position", "Z-position", "Both positions"), table_storage)
write.table(table_storage, file = "STable2_H_Q2.csv", sep = ",", quote = F, row.names = F, col.names = F)

################################################################################
# Make Q2 plot for all features (Figure 5)
xlim = c(0, 0.8)
ylim = c(0, 0.8)
marker_type = "o"
marker_type2 = 18

pdf("figure_dump/Fig5_Predicting_SN.pdf")
### x position
plot(1, type = "n", pch = marker_type, col = "red", xlab = "Kac", ylab = "K", xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")

colors = c("red", "blue", "purple")
for (n_position in 1:3) {
  
  par(new = T)
  # plot individuals
  plot(Q2_int[[1,2]][1:39 + (n_position - 1) * 39], Q2_int[[2,2]][1:39 + (n_position - 1) * 39], pch = marker_type, col = colors[[n_position]], xlab = "", ylab = "", xlim = xlim, ylim = ylim, axes = F, cex = 0.75, xaxt = "n", yaxt = "n")
  
  par(new = T)
  # plot all properties
  plot(Q2_int[[1,1]][[n_position]], Q2_int[[2,1]][[n_position]], pch = marker_type2, col = colors[[n_position]], xlab = "", ylab = "", xlim = xlim, ylim = ylim, axes = F, cex = 3, xaxt = "n", yaxt = "n")
 
}
axis(1, at = seq(0, 0.8, 0.2), labels = seq(0, 0.8, 0.2))
axis(2, at = seq(0, 0.8, 0.2), labels = seq(0, 0.8, 0.2))
lines(x = xlim, y = ylim, lty = 3)

# plot random (choose all features and both positions)
par(new = T)
plot(Q2_int_random[[1,1]][[3]], Q2_int_random[[2,1]][[3]], pch = marker_type2, col = "grey", xlab = "", ylab = "", xlim = xlim, ylim = ylim, axes = F, cex = 3, xaxt = "n", yaxt = "n")

dev.off()



