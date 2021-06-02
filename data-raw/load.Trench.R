## load.Trench.R: Process the 2015 trench oxygen isotope data
##
## Author: Thomas Münch (thomas.muench@awi.de), Alfred-Wegener-Institut, 2018
##

library(usethis)

# 'TrenchR' package to load the data
library(TrenchR)

# Approximate data onto even vertical resolution of 3 cm
res <- 3
depth <- seq(min(t15.trench1$depth), max(t15.trench1$depth), res)

t15.1 <- apply(t15.trench1$data[1, , ], 2, function(prf) {
    approx(t15.trench1$depth, prf, depth)$y})

t15.2 <- apply(t15.trench2$data[1, , ], 2, function(prf) {
    approx(t15.trench1$depth, prf, depth)$y})

# Cut upper incomplete region
t15.2 <- t15.2[-(1 : 6), ]
t15.1 <- t15.1[-(1 : 6), ]

# Combine all trench profiles into one data frame
t15 <- as.data.frame(cbind(t15.1, t15.2))
# Document time of processing
attr(t15, "info") <- sprintf("Processed on: %s.", Sys.time())

usethis::use_data(t15, overwrite = TRUE)

