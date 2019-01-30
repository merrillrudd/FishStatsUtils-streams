
#' @title
#' Check predicted encounter probability against observed encounter frequency
#'
#' @description
#' \code{plot_encounter_diagnostic} is a diagnostic function for checking validity of the encounter-probability component of a spatio-temporal model
#'
#' @inheritParams plot_maps
#' @inheritParams plot_biomass_index
#' @param ... arguments passed to \code{par}
#' @importFrom FishStatsUtils plot_lines
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{Diag_i}{Diagnostic output for each sample \code{i}}
#'   \item{Diag_z}{Diagnostic output for each bin \code{z}}
#' }

#' @export
plot_encounter_diagnostic = function( Report, Data, cutpoints_z=seq(0,1,length=21), interval_width=1.96, savedir=paste0(getwd(),"/"),
  PlotName="Diag--Encounter_prob.png", ... ){

  # Get bin for each datum
  ## Report$R1_i = probability of occurrence
  ## place the predicted probability of occurrence for each observation i into bins in 0.05 intervals between [0,1]
  z_i = cut( Report$R1_i, breaks=cutpoints_z, include.lowest=TRUE )
  midpoints_z = rowMeans( cbind(cutpoints_z[-1],cutpoints_z[-length(cutpoints_z)]) )

  # Get encounter frequency for each bin
  ## each observation had either encounter or non-encounter
  freq_z = tapply( ifelse(Data$b_i>0,1,0), INDEX=z_i, FUN=mean )

  # Get expectation given model
  num_z = tapply( Report$R1_i, INDEX=z_i, FUN=length )
  mean_z = tapply( Report$R1_i, INDEX=z_i, FUN=mean )
  var_z = tapply( Report$R1_i, INDEX=z_i, FUN=function(vec){sum(vec*(1-vec))} )
  sd_mean_z = sqrt(var_z / num_z^2)

  # Plot
  Par = list( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02, ... )
  if(is.null(savedir)==FALSE) png( file=paste0(savedir,"/",PlotName), width=5, height=5, res=200, units="in")
    par( Par )
    plot( x=midpoints_z, y=freq_z, pch=20, cex=1.2, xlim=c(0,1), ylim=c(0,1), xlab="Predicted encounter probability", ylab="Observed encounter frequency" )
    plot_lines( x=midpoints_z[which(!is.na(mean_z))], y=mean_z[which(!is.na(mean_z))], ybounds=(mean_z%o%c(1,1)+sd_mean_z%o%c(-interval_width,interval_width))[which(!is.na(mean_z)),], lwd=2, bounds_type="shading", col_bounds=rgb(1,0,0,0.2), col="red" )
    abline(a=0, b=1, lty="dotted", lwd=2 )
    legend( "topleft", legend=c("Observed","Predicted"), fill=c("black","red"), bty="n")
  if(is.null(savedir)==FALSE) dev.off()

    # df <- data.frame('bin'=midpoints_z, 
    #                 'observed_frequency'=freq_z, 
    #                 "predicted_frequency"=mean_z, 
    #                 'predicted_lcl'=mean_z - interval_width*sd_mean_z, 
    #                 "predicted_ucl"=mean_z + interval_width*sd_mean_z) %>%
    #       na.omit()
    # p <- ggplot(df) +
    #   geom_line(aes(x = bin, y = predicted_frequency), col='red', lwd=2) +
    #   geom_ribbon(aes(bin, ymin = predicted_lcl, ymax = predicted_ucl), fill="red", alpha=0.3) +
    #   geom_point(aes(x = bin, y = observed_frequency), cex=2) +
    #   geom_abline(aes(slope = 1, intercept = 0), linetype=2, lwd=1) +
    #   xlab("Predicted encounter probability") + ylab("Observed encounter probability") +
    #   coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    #   mytheme()

  # Return stuff
  Return = NULL
  Return[["Diag_i"]] = cbind("b_i"=Data$b_i, "z_i"=z_i)
  Return[["Diag_z"]] = cbind("midpoints_z"=midpoints_z, "freq_z"=freq_z, "num_z"=num_z, "mean_z"=mean_z, "var_z"=var_z, "sd_mean_z"=sd_mean_z)
  return( invisible(Return) )
}

