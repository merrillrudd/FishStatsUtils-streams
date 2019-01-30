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
category_names <- unique(obs$survey)

##################################
## model - density observations
##################################

##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = obs$density, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Category" = obs$surveynum)


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

#####################################
## spatial variation
#####################################
FieldConfig = c("Omega1"="IID", "Epsilon1"=0, "Omega2"="IID", "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,1)

Data = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
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

# Change starting values for kappa - was leading to decorrelation at much smaller distances than the minimum distance between sites
Obj$par[grep("logkappa",names(Obj$par))]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

# check that all parameters have a reasonable gradient
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Initial run and hessian check
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=NULL, bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
H = optimHess( par=Opt1$par, fn=Obj$fn, gr=Obj$gr )
check <- TMBhelper::Check_Identifiable(Obj)

# Re-run from last MLE
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=NULL, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]


## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data, savedir=NULL)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, savedir=NULL, plot=2) 

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=2 )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

StreamUtils::plot_maps(plot_set=3, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, cex=0.5, alpha=0.8)

StreamUtils::plot_maps(plot_set=12, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, cex=0.5, alpha=0.8)

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, savedir=NULL, strata_names="Coho salmon", "category_names"=category_names )

FishStatsUtils::plot_range_index(Sdreport=Opt[["SD"]], Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=NULL)


#####################################
## spatiotemporal variation
#####################################

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
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

# Change starting values for kappa - was leading to decorrelation at much smaller distances than the minimum distance between sites
Obj$par[grep("logkappa",names(Obj$par))]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

# check that all parameters have a reasonable gradient
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Initial run and hessian check
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=NULL, bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
ParHat1 = Obj$env$parList()
#H = optimHess( par=Opt1$par, fn=Obj$fn, gr=Obj$gr )

################
# I see that L_epsilon1_z[2] and L_epsilon2_z[1] and L_epsilon2_z[2] are going to zero, so I turn them off
################
FieldConfig[4] = 0
RhoConfig[4] = 0
Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1,
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
Map = TmbList$Map
Map$L_epsilon1_z = factor( c(1,NA) )
Params = TmbList$Parameters
Params$L_epsilon1_z[2] = 0
TmbList = Build_TMB_Fn("TmbData"=Data, "Parameters"=Params, "Version"=Version, "Map"=Map, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

# Re-run optimizer
Opt2 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=NULL, bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
ParHat2 = Obj$env$parList()
H = optimHess( par=Opt2$par, fn=Obj$fn, gr=Obj$gr )

# Re-run from last MLE
## may be sensitive to starting values***
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt2$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=NULL, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

Report = Obj$report()

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data, savedir=NULL)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, savedir=NULL, plot=2) 

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=2 )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

StreamUtils::plot_maps(plot_set=3, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, cex=0.5, alpha=0.8)

StreamUtils::plot_maps(plot_set=12, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, cex=0.5, alpha=0.8)

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, savedir=NULL, strata_names="Coho salmon", "category_names"=category_names )
