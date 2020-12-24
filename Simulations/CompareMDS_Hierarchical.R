# install.packages("readxl")
# install.packages("dendextend")
library(readxl)
library(dendextend)
#library("readxl")
setwd("/Users/annalittle/MDS/Simulations")
k5_n200_sigp2 <- read_excel("k5_n200_sigp2.xlsx", col_names = FALSE)
View(k5_n200_sigp2)
Labels <- c(rep(1,200), rep(2,200), rep(3,200), rep(4,200), rep(5,200))
par(mfrow = c(1,2))
# Colors based on 5 clusters from the tree:
dend <- as.dendrogram(hclust(dist(k5_n200_sigp2),method = "complete"))
dend <- color_labels(dend, k = 5)
plot(dend)
title(main = "Colors based on Hierarchical Clustering")
# Colors based on the labels:
dend <- as.dendrogram(hclust(dist(k5_n200_sigp2),method = "complete"))
ReorderedLabels <- Labels[order.dendrogram(dend)]
labels_colors(dend) <- ReorderedLabels
plot(dend)
title(main = "Colors based on Labels")

k5_n200_sigp3 <- read_excel("k5_n200_sig.3.xlsx", col_names = FALSE)
View(k5_n200_sigp3)
Labels <- c(rep(1,200), rep(2,200), rep(3,200), rep(4,200), rep(5,200))
par(mfrow = c(1,2))
# Colors based on 5 clusters from the tree:
dend <- as.dendrogram(hclust(dist(k5_n200_sigp3),method = "complete"))
dend <- color_labels(dend, k = 5)
plot(dend)
title(main = "Colors based on Hierarchical Clustering")
# Colors based on the labels:
dend <- as.dendrogram(hclust(dist(k5_n200_sigp3),method = "complete"))
ReorderedLabels <- Labels[order.dendrogram(dend)]
labels_colors(dend) <- ReorderedLabels
plot(dend)
title(main = "Colors based on Labels")