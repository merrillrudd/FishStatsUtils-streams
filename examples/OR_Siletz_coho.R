rm(list=ls())

#################
## Packages
#################

devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("merrillrudd/RuddR")
devtools::install_github("merrillrudd/StreamUtils")

library(VAST)
library(StreamUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

#########################
## read in data
#########################

data <- data("or_siletz_coho", package="StreamUtils")

network <- or_siletz_coho[["network"]]
obs <- or_siletz_coho[["observations"]]

Network_sz = network %>% select( -c("long","lat") )
