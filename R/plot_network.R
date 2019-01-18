#' Plot stream network and optional observations
#' 
#' \code{plot_network} generates figure showing network nodes, directional arrows, and observations
#' 
#' @param Extrapolation_List Output from \code{Prepare_Extrapolation_Data_Fn}
#' @param Spatial_List Output from \code{FishStatsUtils::make_spatial_list}
#' @param Data output of VAST::Data_Fn which includes parent_s and child_s
#' @param Data_Geostat data-frame of data (with columns 'Year', 'Lon', 'Lat' or 'E_km', 'N_km' at a minimum)
#' @param PlotDir Directory for plots
#' @param observations default TRUE to include observations on the figure, FALSE to show only network
#' @param arrows default TRUE to show directional connections between nodes, but FALSE may prevent crash if there are too many points. 
#' @param root default TRUE to show root nodes, FALSE in case there are other root nodes that are not meaningful.
#' @param ... addition inputs to \code{plot}
#' 
#' @return Figure plotting stream network and observations

#' @export

plot_network <- function(Extrapolation_List, Spatial_List, Data, Data_Geostat, observations=TRUE, arrows=TRUE, root=TRUE, PlotDir=NULL, ...){

	## observation locations
	obs1 <- data.frame(Extrapolation_List$Data_Extrap)
	if(all(Data_Geostat$Lat==obs1$Lat) & all(Data_Geostat$Lon==obs1$Lon)){
		obs <- cbind.data.frame(Data_Geostat, "E_km"=obs1$E_km, "N_km"=obs1$N_km)
	} else {
		if(grepl("E_km", names(Data_Geostat)) & grepl("N_km", names(Data_Geostat))){
			obs <- Data_Geostat
		} else{
			stop("Extrapolation_List and Data_Geostat latitude and longitude locations do not match -- this function requires both")
		}
	}

	Network <- cbind.data.frame('parent_s'=Data$parent_s+1, 'child_s'=Data$child_s+1, Spatial_List$loc_x)

	years <- unique(Data_Geostat$Year)[order(unique(Data_Geostat$Year))]
	Network <- lapply(1:length(years), function(x){
		out <- cbind.data.frame(Network, "Year"=years[x])
		return(out)
	})
	Network <- do.call(rbind, Network)

	### across years
	aa <- ggplot(Network) +
			mytheme()

	## add roots underneath points
	if(root == TRUE) aa <- aa + geom_point(data = Network %>% filter(parent_s==0), aes(x = E_km, y = N_km), color="goldenrod", cex=5)

	## add points and complete figure
	aa <- aa +
	    geom_point(data = Network %>% filter(Year==years[1]), aes(x = E_km, y = N_km), color = "black", alpha=0.6) +
	    xlab("Eastings") + ylab("Northings")

	## option to add arrows
	if(arrows == TRUE){
		l2 <- lapply(1:nrow(Network), function(x){
			parent <- Network$parent_s[x]
			find <- Network %>% filter(child_s == parent)
			if(nrow(find)>0) out <- cbind.data.frame(Network[x,], 'E2'=find$E_km, 'N2'=find$N_km)
			if(nrow(find)==0) out <- cbind.data.frame(Network[x,], 'E2'=NA, 'N2'=NA)
			# if(nrow(find)>0) out <- cbind.data.frame(Network[x,], 'long2'=find$long, 'lat2'=find$lat)
			# if(nrow(find)==0) out <- cbind.data.frame(Network[x,], 'long2'=NA, 'lat2'=NA)
			return(out)
		})
		l2 <- do.call(rbind, l2)
		aa <- aa + geom_segment(data=l2, aes(x = E2,y = N2, xend = E_km, yend = N_km), arrow=arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), col="gray", alpha=0.6)
	}

	## option to add observations
	if(observations == TRUE) aa <- aa + geom_point(data = obs, aes(x = E_km, y = N_km), color=brewer.pal(3, "Set1")[1], cex=3)

	### by year
	bb <- ggplot(Network) +
			facet_wrap(~Year) +
			mytheme() 

	## add roots underneath points
	if(root == TRUE) bb <- bb + geom_point(data = Network, aes(x = E_km, y = N_km), color="goldenrod", cex=5)

	## add points and complete figure
	bb <- bb +
	    geom_point(data = Network, aes(x = E_km, y = N_km), color = "black", alpha=0.6) +
	    xlab("Eastings") + ylab("Northings")

	## option to add arrows
	if(arrows == TRUE){
		l2 <- lapply(1:nrow(Network), function(x){
			parent <- Network$parent_s[x]
			find <- Network %>% filter(child_s == parent)
			if(nrow(find)>0) out <- cbind.data.frame(Network[x,], 'E2'=find$E_km, 'N2'=find$N_km)
			if(nrow(find)==0) out <- cbind.data.frame(Network[x,], 'E2'=NA, 'N2'=NA)
			# if(nrow(find)>0) out <- cbind.data.frame(Network[x,], 'long2'=find$long, 'lat2'=find$lat)
			# if(nrow(find)==0) out <- cbind.data.frame(Network[x,], 'long2'=NA, 'lat2'=NA)
			return(out)
		})
		l2 <- do.call(rbind, l2)
		bb <- bb + geom_segment(data=l2, aes(x = E2,y = N2, xend = E_km, yend = N_km), arrow=arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), col="gray", alpha=0.6)
	}

	## option to add observations
	if(observations == TRUE) bb <- bb + geom_point(data = obs, aes(x = E_km, y = N_km), color=brewer.pal(3, "Set1")[1], cex=3)

	if(all(is.null(PlotDir))==FALSE){
		ggsave(file.path(PlotDir, "Network.png"), aa)
		ggsave(file.path(PlotDir, "Network_byYear.png"), aa)
	}
	print(aa)
	dev.new()
	print(bb)
}