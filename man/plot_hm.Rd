% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plot_hm}
\alias{plot_hm}
\title{A Function That Plots Heatmaps of Differentialy Expressed Genes Accros Conditions}
\usage{
plot_hm(exp_grp_key, case_fctr, anova_mw_res, study, ctrl_fctr, main,
  mw_pval = "mw_pval_adj", gene_pf_colname = "gene_name", key,
  PLOT_GS = FALSE, fc_thres = c(-1.5, -1.2, 1.2, 1.5),
  PLOT_PVAL = TRUE)
}
\arguments{
\item{exp_grp_key}{the column name of the exprimental grouping to considere}

\item{case_fctr}{The label of the group that we considere as differentialy expressed}

\item{anova_mw_res}{A dataframe with ANOVA and Mann-Whitney results}

\item{study}{An Reference classes object that contains data and metadata}

\item{ctrl_fctr}{The label of the group that we considere as the reference}

\item{main}{A character string to explicit the title of the plot}

\item{mw_pval}{the the suffix the column name corresponding to the Mann-Whitney p-value used to}

\item{gene_pf_colname}{A character string specifying the name of the platform column to use for the gene name}

\item{key}{A character string that will suffix the column names of the resulting data frame}

\item{PLOT_GS}{A boolean set to TRUE if genes names to be ploted}

\item{fc_thres}{A vector of integer coprresponding to the foldchange thresholds}

\item{PLOT_PVAL}{A boolean set to TRUE if ANOVA and Mann-Whitney need to be ploted}
}
\value{
probes used to plot the heatmap
}
\description{
This function plots heatmaps of differentialy expressed genes across conditions.
}
