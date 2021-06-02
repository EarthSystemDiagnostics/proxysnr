## clim-site-par.R: Read in site parameters to provide in package
##
## Author: Thomas Münch (thomas.muench@awi.de), Alfred-Wegener-Institut, 2018
##

library(usethis)

dml <- read.table("data-raw/dml_diffusion-par.txt",
                  header = TRUE, skip = 10, sep = "\t")
wais <- read.table("data-raw/wais_diffusion-par.txt",
                   header = TRUE, skip = 14, sep = "\t")

clim.site.par <- rbind(dml, wais)

usethis::use_data(clim.site.par, overwrite = TRUE)

