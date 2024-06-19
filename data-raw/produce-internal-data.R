## Produce internal proxysnr data:
##
## - clim.site.par: climatological parameters for diffusion length calculations
## - diffusion.length: calculated diffusion lengths for the DML and WAIS sites
## - t15.noise: DML 2015 trench oxygen isotope noise spectra
##
## Author: Thomas Münch (thomas.muench@awi.de),
##         Alfred-Wegener-Institut, 2018 (update 2024)
##

# set working directory to your package source folder
# setwd("<path>")

# require package 'FirnR' for density and diffusion length calculations;
# available on request from the `proxysnr` authors
library(FirnR)

# require package 'TrenchR' for access to t15 trench data
# remotes::install_github("EarthSystemDiagnostics/TrenchR")
library(TrenchR)

library(usethis)

# ==============================================================================
# I. Calculate DML and WAIS diffusion lengths
# ==============================================================================

#' Local climatological parameters for diffusion length calculations
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{core:}{the name of the firn-core site}
#'   \item{temperature:}{10-m firn temperature (or annual mean air
#'     temperature), in degree Celsius}
#'   \item{acc.rate:}{accumulation rate, in kg/m^2/yr}
#'   \item{pressure:}{surface air pressure, in mbar}
#'   \item{height:}{surface elevation, in m}
#' }
#'
#' @source
#' All values for firn cores in the rows 1-15 from Oerter et al. (2000) (Table
#' 1, 4), except: temperature for FB9807 adapted from B32 (1 km), temperature
#' for FB9815 adapted from FB9805 (28 km), temperature for FB9811 adapted from
#' FB9812 (19 km).
#'
#' For firn cores in the rows 16-20, accumulation rates and elevation data are
#' for core WDC2005A taken from WAIS Divide Project Members (2013) (p.440) and
#' for the other cores from Kaspari et al. (2004) (Table 1), respectively.
#' WDC2005A temperature is taken from WAIS Divide Project Members (2013)
#' (p.440), the other temperature values are based on ERA-Interim mean
#' temperature anomalies (Dee et al., 2011) with respect to WDC2005A.
#'
#' For all cores, the pressure values are calculated with the barometric formula
#' using the local elevation and temperature data.
#'
#' @references
#' Dee, D. P., et al.: The ERA-Interim reanalysis: configuration and performance
#' of the data assimilation system, Q. J. R. Meteorol. Soc., 137(656), 553–597,
#' \url{https://doi.org/10.1002/qj.828}, 2011.
#' 
#' Kaspari, S., et al.: Climate variability in West Antarctica derived from
#' annual accumulation-rate records from ITASE firn/ice cores, Ann. Glaciol.,
#' 39(1), 585-594, \url{https://doi.org/10.3189/172756404781814447}, 2004.
#'
#' Münch, T. and Laepple, T.: What climate signal is contained in
#' decadal- to centennial-scale isotope variations from Antarctic ice cores?
#' Clim. Past, 14, 2053–2070, \url{https://doi.org/10.5194/cp-14-2053-2018},
#' 2018.
#'
#' Oerter, H., et al.: Accumulation rates in Dronning Maud Land, Antarctica, as
#' revealed by dielectric-profiling measurements of shallow firn cores,
#' Ann. Glaciol., 30(1), 27-34, 2000.
#'
#' WAIS Divide Project Members: Onset of deglacial warming in West Antarctica
#' driven by local orbital forcing, Nature, 500(7463), 440–444,
#' \url{https://doi.org/10.1038/nature12376}, 2013.
#'
#' @author Thomas Münch
#'
readClimParameters <- function() {

  dml <- read.table("data-raw/dml_diffusion-par.txt",
                    header = TRUE, skip = 10, sep = "\t")
  wais <- read.table("data-raw/wais_diffusion-par.txt",
                     header = TRUE, skip = 14, sep = "\t")

  rbind(dml, wais)

}

# Parameters to simulate diffusion length

# climatological site parameters
clim.site.par <- readClimParameters()

# simulated core length [m]
core.len <- 500
# vertical resolution [m]
z.res <- 1/100
# depth profile
depth <- seq(z.res, core.len, z.res)
# time resolution [yr]
t.res <- 1/2
# surface firn density [kg/m^3]
rho.surface <- 340

# structure to store results
diffusion.length <- list()

#-------------------------------------------------------------------------------
# Simulation for DML1 cores

# Relevant site parameters
param <- clim.site.par[1 : 15, ]

# Number of cores
nc <- length(dml$dml1)

# Number of simulated ages
nt <- length(dml$dml1[[1]]) / t.res

# Array to store diffusion length in years
sigma <- array(dim = c(nt, nc + 1))

# Loop over individual core sites
for (i in 1 : nc) {

  # Herron-Langway density profile for site 'i'
  HL  <- FirnR::DensityHL(depth = depth,
                          rho.surface = rho.surface,
                          T = param$temperature[i] + 273.15,
                          bdot = param$acc.rate[i])

  # age of core according to Herron-Langway solution and const. acc. rate
  t <- HL$depth.we * 1000 / param$acc.rate[i]

  # diffusion length in [cm] for site 'i' as a function of depth
  sig.cm <- FirnR::DiffusionLength(depth = depth,
                                   rho = HL$rho,
                                   T = param$temperature[i] + 273.15,
                                   P = param$pressure[i],
                                   bdot = param$acc.rate[i])

  # convert diffusion length from [cm] to [yr]
  sig.yr <- 1e-2 * sig.cm * (HL$rho / param$acc.rate[i])
  
  # diffusion length in [yr] on equidistant time grid
  sigma[, 1]     <- seq_len(nt) * t.res
  sigma[, i + 1] <- approx(t / t.res, sig.yr, seq_len(nt), rule = 2)$y
  
}

# Store results
sigma <- as.data.frame(sigma)
colnames(sigma) = c("Time", as.character(param$core))
diffusion.length$dml1 <- sigma


#-------------------------------------------------------------------------------
# Simulation for DML2 cores

param <- clim.site.par[1 : 3, ]
nc <- length(dml$dml2)
nt <- length(dml$dml2[[1]]) / t.res

sigma <- array(dim = c(nt, nc + 1))

for (i in 1 : nc) {

  HL  <- FirnR::DensityHL(depth = depth,
                          rho.surface = rho.surface,
                          T = param$temperature[i] + 273.15,
                          bdot = param$acc.rate[i])

  t <- HL$depth.we * 1000 / param$acc.rate[i]

  sig.cm <- FirnR::DiffusionLength(depth = depth,
                                   rho = HL$rho,
                                   T = param$temperature[i] + 273.15,
                                   P = param$pressure[i],
                                   bdot = param$acc.rate[i])

  sig.yr <- 1e-2 * sig.cm * (HL$rho / param$acc.rate[i])
  
  sigma[, 1]     <- seq_len(nt) * t.res
  sigma[, i + 1] <- approx(t / t.res, sig.yr, seq_len(nt), rule = 2)$y
  
}

sigma <- as.data.frame(sigma)
colnames(sigma) = c("Time", as.character(param$core))
diffusion.length$dml2 <- sigma


#-------------------------------------------------------------------------------
# Simulation for WAIS cores

param <- clim.site.par[16 : 20, ]
nc <- length(wais)
nt <- length(wais[[1]]) / t.res

sigma <- array(dim = c(nt, nc + 1))

for (i in 1 : nc) {

  HL  <- FirnR::DensityHL(depth = depth,
                          rho.surface = rho.surface,
                          T = param$temperature[i] + 273.15,
                          bdot = param$acc.rate[i])

  t <- HL$depth.we * 1000 / param$acc.rate[i]

  sig.cm <- FirnR::DiffusionLength(depth = depth,
                                   rho = HL$rho,
                                   T = param$temperature[i] + 273.15,
                                   P = param$pressure[i],
                                   bdot = param$acc.rate[i])

  sig.yr <- 1e-2 * sig.cm * (HL$rho / param$acc.rate[i])
  
  sigma[, 1]     <- seq_len(nt) * t.res
  sigma[, i + 1] <- approx(t / t.res, sig.yr, seq_len(nt), rule = 2)$y
  
}

sigma <- as.data.frame(sigma)
colnames(sigma) = c("Time", as.character(param$core))
diffusion.length$wais <- sigma

# ==============================================================================
# II. Calculate DML 2015 trench oxygen isotope noise spectra
# ==============================================================================

#' Calculate trench noise spectrum given a constant accumulation rate
#'
#' @param acc.rate assumed mean accumulation rate at the trench site [cm/yr] to
#'   convert the depth axis into a time axis.
#' @param res depth resolution of the trench data; defaults to 3 cm.
#' @param neff effective number of trench records to account for the spatial
#'   autocorrelation of the stratigraphic noise; defaults to 19 for the total
#'   available number of 22 profiles (see Münch et al., 2017,
#'   https://doi.org/10.5194/tc-11-2175-2017).
#' @param df.log Gaussian kernel width in log space to smooth the spectral
#'   estimate; defaults to 0.1.
#' @return a list of class \code{"spec"} with two elements: the frequency axis
#'   (element `freq`) and the estimate noise PSD (element `spec`).
#'
#' @author Thomas Münch
#'
getTrenchNoise <- function(t15, acc.rate = 25, res = 3,
                           neff = 19, df.log = 0.1) {

  SeparateSignalFromNoise(
    ObtainArraySpectra(
      t15, res = res / acc.rate, neff = neff, df.log = df.log))$noise

}

# approximate trench data onto even vertical resolution of 3 cm
res <- 3
depth.orig <- TrenchR::getZ(TrenchR::t15.trench1) # same as for trench2
depth <- seq(min(depth.orig), max(depth.orig), res)

t15.1 <- apply(TrenchR::make2D(TrenchR::t15.trench1, simplify = TRUE), 2,
               function(prf) {
                 approx(depth.orig, prf, depth)$y
               })
t15.2 <- apply(TrenchR::make2D(TrenchR::t15.trench2, simplify = TRUE), 2,
               function(prf) {
                 approx(depth.orig, prf, depth)$y
               })

# cut upper incomplete region and neglect excess T15-1 profile #7
t15.1 <- t15.1[-(1 : 6), -7]
t15.2 <- t15.2[-(1 : 6), ]

# combine all trench profiles into one data frame
t15 <- as.data.frame(cbind(t15.1, t15.2))

# estimate noise spectra for three different accumulation rates

t15.noise <- list(

  # lower bound accumulation rate
  lower = getTrenchNoise(t15, acc.rate = 25 - 5),
  # mean accumulation rate
  mean = getTrenchNoise(t15),
  # upper bound accumulation rate
  upper = getTrenchNoise(t15, acc.rate = 25 + 5)

)

# ==============================================================================
# Save package data
# ==============================================================================

usethis::use_data(clim.site.par, diffusion.length, t15.noise,
                  internal = TRUE, overwrite = TRUE)
