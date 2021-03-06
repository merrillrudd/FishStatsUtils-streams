% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_quantile_diagnostic.R
\name{plot_quantile_diagnostic}
\alias{plot_quantile_diagnostic}
\title{Diagnostic QQ function}
\usage{
plot_quantile_diagnostic(TmbData, Report, DateFile = paste0(getwd(),
  "/"), savedir = paste0(DateFile, "/QQ_Fn/"),
  FileName_PP = "Posterior_Predictive",
  FileName_Phist = "Posterior_Predictive-Histogram",
  FileName_QQ = "Q-Q_plot", FileName_Qhist = "Q-Q_hist", plot = 1:4,
  category_names = NULL)
}
\arguments{
\item{TmbData}{TMB Model input data list}

\item{Report}{TMB Model output data list}

\item{savedir}{Directory to save plots, if NULL is specified then do not save plot}

\item{FileName_PP}{If NULL is specified then do not save this type of plot}

\item{FileName_Phist}{If NULL is specified then do not save this type of plot}

\item{FileName_QQ}{If NULL is specified then do not save this type of plot}

\item{FileName_Qhist}{If NULL is specified then do not save this type of plot}

\item{plot}{default 1:4 to plot all diagnostics}

\item{category_names}{names of categories for plotting labels}
}
\value{
A list containing results for each specified categories
}
\description{
Diagnostic QQ function
}
\examples{
Q <- QQ_Fn(TmbData = TmbData, Report = Report)
}
