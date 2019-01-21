rm(list=ls())

#################
## Packages
#################

devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("merrillrudd/FishStatsUtils")
devtools::install_github("merrillrudd/RuddR")
devtools::install_github("merrillrudd/StreamUtils")

library(VAST)
library(FishStatsUtils)
library(StreamUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

#########################
## read in data
#########################

data <- data("nz_waitaki_longfin_eel", package="StreamUtils")

network <- nz_waitaki_longfin_eel[["network"]]
obs <- nz_waitaki_longfin_eel[["observations"]]

Network_sz = network %>% select( -c("long","lat") )

##################################
## model - encounter observations
##################################

##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0)


## include latitude and longitude for user-supplied area
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs$lat, 
                    "Lon"=obs$long,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$lat, 
                          "LON_intensity"=network$long, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

########################
## No spatial variation
########################
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, Data=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=FALSE, root=FALSE)

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data, DirName=NULL)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, save_dir=NULL, plot=2) #StreamUtils::

## Plot Pearson residuals
plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

##### Model output
## density surface for each year
Dens_xt <- plot_maps(plot_set=3, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", savedir=NULL)

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, DirName=NULL, strata_names = "Longfin eels" )

##################
# Now try turning on spatial variation
##################
FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)

OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0,
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=FALSE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data, DirName=NULL)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, save_dir=NULL, plot=2) #StreamUtils::

## Plot Pearson residuals
plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

##### Model output
## density surface for each year
Dens_xt <- plot_maps(plot_set=3, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", savedir=NULL, cex=0.5)

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, DirName=NULL, strata_names = "Longfin eels" )
