## calculate-diffusion-length.R:
## Calculate diffusion lengths for the DML and WAIS firn-core sites
##
## Author: Thomas Münch (thomas.muench@awi.de), Alfred-Wegener-Institut, 2018
##

library(usethis)

# Package needed for climatological site and core parameters 
library(proxysnr)

# Package 'FirnR' for density and diffusion length calculations;
# available on request from the authors
library(FirnR)

# Parameters to simulate diffusion length

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

# Structure to store results
diffusion.length <- list()


#-------------------------------------------------------------------------------
# Simulation for DML1 cores

# Relevant site parameters
param <- proxysnr::clim.site.par[1 : 15, ]

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

param <- proxysnr::clim.site.par[1 : 3, ]
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

param <- proxysnr::clim.site.par[16 : 20, ]
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


#-------------------------------------------------------------------------------

usethis::use_data(diffusion.length, overwrite = TRUE)


