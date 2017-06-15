
rm(list=ls())

setwd("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/")

source("~/Dropbox/Albert Xue/Manuscripts/ACS Analytical Chemistry 2016/Code/support_functions/")

# load experimental data from Kuo et al. 
load("../Data/Kuo_KDAC_experimental_data.RData")
# indices that correspond to different types of replicates
# there is a pair of biological replicates for each condition/timepoint
indices = list(list(1:2, 3:4), list(5:6, 7:8), list(9, 10))

# load signal to noise values
load("../Data/K_array_data.RData")
SN_K = colMeans(SN_K)

##############################################################################
# plotting

xlim = c(0, 250)
ylim = c(0, 0.18)
xlab = "Mean signal fidelity of quantiles"
ylab = "Standard deviation"
pch_list = rep(c(0,2), each = 3)
color_list = rep(c("blue", "green", "red"), 2)

# storage for correlations
corr_result = matrix(NA, nrow = 2, ncol = 3)
corr_result_random = matrix(NA, nrow = 2, ncol = 3)

pdf(file = "figure_dump/SFig3_SN_correlation.pdf")

# indices for the five quintiles
quintile_indices = split(1:361, cut(seq_along(1:361), 5, labels = FALSE))
# SN values of the five quintiles
quintile_SN = sapply(quintile_indices, function(i) mean(sort(SN_K, decreasing = T)[i]))

plot(1, type = "n", xlim = xlim, ylim = ylim, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
par(new = T)
counter = 0
for (n_day in 1:2) {
  for (n_cond in 1:3) {
    counter = counter + 1
    
    # grab data for specified day/condition
    day_data = matrix(data_mat[indices[[n_cond]][[n_day]],], ncol = 361)
    # sort each by SN
    day_data = matrix(day_data[,order(SN_K, decreasing = T)], ncol = 361)
    
    # standard deviation of peptides in the five quintiles
    y_values = sapply(1:5, function(i) sd(day_data[, quintile_indices[[i]]]))
    
    plot(x = quintile_SN, y = y_values, xlim = xlim, ylim = ylim, axes = F, xlab = "", ylab = "", pch = pch_list[[counter]], col = color_list[[counter]], cex = 1.5, bty = "n", lwd = 2)
    par(new = T)
    
    # correlation between peptide sd and mean SN of each quantile
    sd_values = sapply(quintile_indices, function(i) sd(day_data[,i]))
    mean_values = sapply(quintile_indices, function(i) mean(sort(SN_K, decreasing=T)[i]))
    corr_result[[n_day, n_cond]] = cor(sd_values, mean_values)
    
    # same for random (do it multiple times)
    temp = c()
    for (n_rand in 1:100) {
      day_data_random = day_data[,sample(1:361), drop = F]
      sd_values = sapply(quintile_indices, function(i) sd(day_data_random[,i]))
      mean_values = sapply(quintile_indices, function(i) mean(sort(SN_K, decreasing=T)[i]))
       temp = c(temp, cor(sd_values, mean_values))
    }
    corr_result_random[[n_day, n_cond]] = mean(temp)
  }
}
axis(side = 1, at = seq(0,200,50), labels=seq(0,200,50))
axis(side = 2, at = seq(0,0.18,0.03), labels=seq(0,0.18,0.03))

legend("bottomleft", legend = c("Day 0", "Day 6"), pch = c(0, 2), cex = 1)
legend("topright", legend = c("Lysate", "TSA", "NIC"), fill = c("blue", "green", "red"))

dev.off()



