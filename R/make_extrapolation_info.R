
#' Build extrapolation grid
#'
#' \code{make_extrapolation_data} builds an object used to determine areas to extrapolation densities to when calculating indices
#'
#' @param strata.limits an input for determining stratification of indices (see example script)
#' @param stream_info data frame with 3 columns: latitude, longitude, and child node of stream segment (upstream node) associated with each observation
#' @inheritParams Convert_LL_to_UTM_Fn

#' @return Tagged list used in other functions
#' \describe{
#'   \item{a_el}{The area associated with each extrapolation grid cell (rows) and strata (columns)}
#'   \item{Data_Extrap}{A data frame describing the extrapolation grid}
#'   \item{zone}{the zone used to convert Lat-Long to UTM by PBSmapping package}
#'   \item{flip_around_dateline}{a boolean stating whether the Lat-Long is flipped around the dateline during conversion to UTM}
#'   \item{Area_km2_x}{the area associated with each row of Data_Extrap, in units square-kilometers}
#' }

#' @export
make_extrapolation_info = function( strata.limits, stream_info=NULL, zone=NA, flip_around_dateline=TRUE, ... ){

  # Infer strata
  if( is.null(strata.limits)){
    strata.limits = data.frame('STRATA'="All_areas")
  }
  message("Using strata ", strata.limits)

  # Read extrapolation data
  Data_Extrap <- stream_info

  # Survey areas
  Area_km2_x = Data_Extrap[,'Area_km2']
  
  # Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_LAT_DD"=Data_Extrap[,'Lat'], "BEST_LON_DD"=Data_Extrap[,'Lon'])
  if( "Depth" %in% colnames(Data_Extrap) ){
    Tmp = cbind( Tmp, Data_Extrap[,'Depth'] )
  }
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits), dimnames=list(NULL,strata.limits[,'STRATA'])))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp, MARGIN=1, FUN=match_strata_fn, strata_dataframe=strata.limits[l,,drop=FALSE])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, Area_km2_x)
  }

  # Convert extrapolation-data to an Eastings-Northings coordinate system
  tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Extrap[,'Lon'], Lat=Data_Extrap[,'Lat'], zone=zone, flip_around_dateline=flip_around_dateline)                                                         #$

  # Extra junk
  Data_Extrap = cbind( Data_Extrap, 'Include'=1)
  if( all(c("E_km","N_km") %in% colnames(Data_Extrap)) ){
    Data_Extrap[,c('E_km','N_km')] = tmpUTM[,c('X','Y')]
  }else{
    Data_Extrap = cbind( Data_Extrap, 'E_km'=tmpUTM[,'X'], 'N_km'=tmpUTM[,'Y'] )
  }

  # Return
  Extrapolation_List = list( "a_el"=a_el, "Data_Extrap"=Data_Extrap, "zone"=attr(tmpUTM,"zone"), "flip_around_dateline"=flip_around_dateline, "Area_km2_x"=Area_km2_x)
  # Return
  return( Extrapolation_List )
}
