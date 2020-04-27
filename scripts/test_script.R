library(ggplot2)
library(gridExtra)
library(nlme)
library(AICcmodavg)
library(car)
library(tidyverse)
library(here)

setwd(here())

############################################
############  mass standardising  ##########
############################################

raw_data <- read.csv("data/raw_metabolic_trait.csv", header = T)
head(raw_data)

## define constants

E <- 0.63  # average activation energy
k <- 8.617333*10^-5  # boltzman constant
Kelvin <- 273.15  # kelvin conversion

## change rates to whole organism metabolic rates
## remove temperature effect from data with Arrhenius function
## calculate the natural logarithm of metabolic rates and mass

raw_data$y.smr <- log((raw_data$smr*raw_data$mass)*exp(E/(k*(raw_data$temp+Kelvin))))
raw_data$y.mmr <- log((raw_data$mmr*raw_data$mass)*exp(E/(k*(raw_data$temp+Kelvin))))
raw_data$ln_mass <- log(raw_data$mass)

## calulate the ln-ln regression of metabolic rates vs mass

smr_scale_mod <- lm(y.smr ~ ln_mass, data = raw_data)
mmr_scale_mod <- lm(y.mmr ~ ln_mass, data = raw_data)

## take the slope of the regressions as scaling coefficients

smr_scale <- as.numeric(smr_scale_mod$coefficients[2])
mmr_scale <- as.numeric(mmr_scale_mod$coefficients[2])

## mass standardise metabolic rates

raw_data$smr_mass_stand <- (raw_data$smr*raw_data$mass)/(raw_data$mass^smr_scale)
raw_data$mmr_mass_stand <- (raw_data$mmr*raw_data$mass)/(raw_data$mass^mmr_scale)


####################################################################
############  model cubic polynomials for metabolic rates ##########
####################################################################

## create a new model dataframe

data <- dplyr::select(raw_data, id, site, temp, smr_mass_stand, mmr_mass_stand)  # select data needed
data$as <- data$mmr_mass_stand - data$smr_mass_stand  # calculate aerobic scope
data$temp0 <- data$temp - 8  # adjust temperature so 8 deg C is the intercept
data <- unite(data, "codetreat", c("site", "temp"), sep = "", remove = FALSE)  # creaet a uniwue code for each temp and location for model weights

## rename columns

data <- rename(data, smr = smr_mass_stand , mmr = mmr_mass_stand)