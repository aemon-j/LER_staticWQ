# Set to your LER folder
setwd('/home/robert/Projects/AEMONJ/LakeEnsemblR')

# Load libraries
library(ggpubr)
library(ggplot2)
library(gotmtools)
library(LakeEnsemblR)
library(tidyverse)
library(pracma)
library(rLakeAnalyzer)
library(LakeMetabolizer)

# your example simulation
setwd('example/feeagh')

ens_out <- 'output/ensemble_output.nc'
vars <- gotmtools::list_vars(ens_out)

p1 <- plot_heatmap(ens_out);p1

# bathymetry
bath <- read_csv('LakeEnsemblR_bathymetry_standard.csv')

# get variable
var_list <- load_var(ncdf = ens_out, var = 'watertemp', return = "list", print = FALSE)
# get depths
deps <- rLakeAnalyzer::get.offsets(var_list[[1]])
# only the selected models
var_list <- var_list[c('GOTM')]
data <- var_list$GOTM[, -c(1)]

# meteorology for wind velocity
start_date <- as.Date(get_yaml_value('LakeEnsemblR.yaml', "time", "start"), ts = 'UTC')
stop_date <- as.Date(get_yaml_value('LakeEnsemblR.yaml', "time", "stop"), ts= 'UTC')
meteo <- read_csv(file = 'LakeEnsemblR_meteo_standard.csv')

id1 <- match(start_date, as.Date(meteo$datetime))
id2 <- match(stop_date, as.Date(meteo$datetime))
meteo.df <- meteo[id1:id2, ]

# get data on even grid
time <- var_list$GOTM[,1]
appr_depths <- seq(from = min(deps), to = max(deps), by = 0.5)
appr_areas <- approx(x = bath$Depth_meter, y = bath$Area_meterSquared, 
                     xout = appr_depths)$y
appr_wtr <- apply(data, 1, function(x) approx(deps, x,appr_depths)$y)
appr_vol <- abs(lead(appr_areas)-appr_areas)
appr_vol[length(appr_vol)] <- appr_areas[length(appr_areas)]
  
# calculate vertical diffusivities
dTdt <- apply(appr_wtr, 1, function(x) abs(lead(x) - x))
dTdt <- t(dTdt)
dTdt[, ncol(dTdt)] <- appr_wtr[, ncol(dTdt)] -appr_wtr[, ncol(dTdt)-1]
dTdz <- apply(appr_wtr, 2, function(x) abs(lead(x) - x))
dTdz[nrow(dTdz), ] <- appr_wtr[nrow(dTdz), ] -appr_wtr[nrow(dTdz)-1, ]
dTdz[which(dTdz == 0)] <- 1e-3

kz_heinz <- matrix(NA, nrow = nrow(appr_wtr), ncol = ncol(appr_wtr))
# Heinz 1990: https://link.springer.com/content/pdf/10.1007/BF00877283.pdf
for (dt in 1: ncol(kz_heinz)) {
  for (dz in 1: nrow(kz_heinz)) {
    kz_heinz[dz, dt] <- trapz(x = appr_depths[dz:nrow(kz_heinz)],
                       dTdt[dz:nrow(kz_heinz), dt] * appr_areas[dz:nrow(kz_heinz)]) /
      (dTdz[dz, dt] * appr_areas[dz])
  }
}

dens <- function(wtemp) {
  density = 999.842594+6.793952e-2*wtemp-9.09529e-3*wtemp^2+1.001685e-4*wtemp^3-1.120083e-6*wtemp^4+6.536332e-9*wtemp^5
  return(density)
}
appr_dens <- dens(appr_wtr)

n = 0.5#0.5
a = 1e-12
n2 <- matrix(NA, nrow = nrow(appr_wtr), ncol = ncol(appr_wtr))
drhodz <- matrix(NA, nrow = nrow(appr_wtr), ncol = ncol(appr_wtr))
for (dt in 1: ncol(n2)) {
  drhodz[, dt] <- lead(appr_dens[, dt]) - appr_dens[, dt]
  drhodz[nrow(drhodz), dt] <- appr_dens[nrow(drhodz), dt] - appr_dens[nrow(drhodz)-1, dt]
  for (dz in 1: nrow(n2)) {
    n2[dz, dt] <- 9.81/appr_dens[dz, dt] * abs(drhodz[dz, dt])
  }
  n2[n2 == 0] <- 1e-10
}
# Quay 1980: http://citeseerx.ist.psu.edu/vie1wdoc/download?doi=10.1.1.568.3595&rep=rep1&type=pdf 
kz_quay <- a * n2^(-n)

# lake number
bth <- as.matrix(bath)
colnames(bth) <- c('depths', 'areas')
ln <- ts.lake.number(wtr = as.data.frame(var_list$GOTM), wnd = data.frame('datetime' = var_list$GOTM$datetime, 
                                                                          'wind' =  
                                                                            meteo.df$Ten_Meter_Elevation_Wind_Speed_meterPerSecond), 
                     wnd.height = 10,
                  bathy = as.data.frame(bth), seasonal = TRUE)
ln$lake.number[is.na(ln$lake.number)] = 1e-1

# Ozkundakci 2011 https://link.springer.com/content/pdf/10.1007/s10750-010-0358-9.pdf
kz_oskundakci <- matrix(NA, nrow = nrow(appr_wtr), ncol = ncol(appr_wtr))
for (dt in 1: ncol(kz_oskundakci)) {
  max_n2 <- max(n2[, dt])
  for (dz in 1: nrow(kz_oskundakci)) {
    kz_oskundakci[dz, dt] <- (200 * n2[dz,dt])/(ln$lake.number[dt]*max_n2) * 1.67*10^(-5)/10000
  }
}

# prevent fluxes to deplete more than what is available
valid_flux <- function(flx, conc){
  if(abs(flx)>abs(conc)){
    if(flx < 0){
      flux = - conc
    }else{
      flux = flx
    }
  }else{
    flux = flx
  }
  return(flux)
}

do <- matrix(NA, nrow = nrow(appr_wtr), ncol = ncol(appr_wtr))
do[, 1] <- 1e4
dx <- mean(diff(appr_depths))
dtt <- as.numeric(mean(diff(time)))

k600 = k600.2.kGAS.base(k600 = k.cole.base(wnd$wind),temperature = appr_wtr[1,], gas = "O2")
o2sat = o2.at.sat.base(temp = appr_wtr[1,], altitude = 300) * 1000

kz = kz_oskundakci 
nep <- 100 #50
sed <- 1000 #1000 
khalf <- 2000
# http://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf 
for (dt in 2: ncol(do)) {
  
  do[1, dt] <- do[1, dt-1] + valid_flux( (1) * kz[1, dt-1]  * dtt / dx**3 * (2 * do[1, dt-1] - 5 * do[1+1, dt-1] + 4 * do[1 +2, dt-1] - do[1+3, dt-1]) +
    k600[dt-1] * (o2sat[dt-1] - do[1, dt-1]) / dx +
    nep * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) -
      # nep  * 1.08^(appr_wtr[dz, dt-1] -20) -
    sed * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) * appr_vol[dz]/(appr_areas[dz] * dx), do[1, dt-1])
  
  do[(nrow(do)), dt] <- do[nrow(do), dt-1] + valid_flux( (1) * kz[nrow(do), dt-1]  * dtt / dx**3 * (2 * do[nrow(do), dt-1] - 5 * do[nrow(do)-1, dt-1] + 4 * do[nrow(do) -2, dt-1] - do[nrow(do)-3, dt-1]) +
    nep * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) -
      # nep  * 1.08^(appr_wtr[dz, dt-1] -20) -
    sed * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) * appr_vol[dz]/(appr_areas[dz] * dx), do[nrow(do), dt-1])
  
  for (dz in 2: (nrow(do)-1)) {
    do[dz, dt] <- do[dz, dt-1] + valid_flux( (1) * kz[dz, dt-1]  * dtt / dx**2 * (do[dz+1, dt-1] - 2 * do[dz, dt-1] + do[dz-1, dt-1]) +
      nep * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) -
        # nep *  1.08^(appr_wtr[dz, dt-1] -20) -
      sed * (do[dz, dt-1]/(do[dz, dt-1] + khalf)) * 1.08^(appr_wtr[dz, dt-1] -20) * appr_vol[dz]/(appr_areas[dz] * dx), do[dz, dt-1])
  }
  do[ which(do[,dt] < 0), dt] <- 0.
}

df.do <- data.frame(cbind(time, t(do)) )
colnames(df.do) <- c("time", as.character(appr_depths))
mdf.do <- reshape2::melt(df.do, "time")
mdf.do$time <- time

df.wtr <- data.frame(cbind(time, data))
colnames(df.wtr) <- c("time", as.character(appr_depths))
mdf.wtr <- reshape2::melt(df.wtr, "time")
mdf.wtr$time <- time

df.kz <- data.frame(cbind(time, t(kz)))
colnames(df.kz) <- c("time", as.character(appr_depths))
mdf.kz <- reshape2::melt(df.kz, "time")
mdf.kz$time <- time

library(RColorBrewer)
library(viridis)
library(patchwork)

p1 <- ggplot(mdf.do, aes(x = time, y = as.numeric(variable), z = value/1000, col = value/1000)) +
  geom_point()  + 
  scale_colour_viridis(option = 'plasma') +
  scale_y_reverse() +
  xlab('time') +
  ylab('depth') +
  labs(colour = 'DO [mg/L]')+
  theme_bw() 

p2 <- ggplot(mdf.wtr, aes(x = time, y = as.numeric(variable), z = value, col = value)) +
  geom_point()  + 
  scale_colour_viridis(option = 'plasma') +
  scale_y_reverse() +
  xlab('time') +
  ylab('depth') +
  labs(colour = 'WTR [degC]')+
  theme_bw()  

p3 <- ggplot(mdf.kz, aes(x = time, y = as.numeric(variable), z = log10(value), col = log10(value))) +
  geom_point()  + 
  scale_colour_viridis(option = 'plasma') +
  scale_y_reverse() +
  xlab('time') +
  ylab('depth') +
  labs(colour = 'log10(Kz) [m2/d]')+
  theme_bw()  


p2 / p1 / p3
