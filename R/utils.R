#' A Function That Plots Survival Curves.
#'
#' This function takes a survival structure such as produced by function Surv of
#' package survival. Takes a (discretized) vector v, indexed as SS.
#' Performs the Cox proportional hazard rate
#' test of SS against v. It plots the survival curves for each
#' level of v. It also prints the p-value of the log-likelihood test.
#'
#' @param SS A survival structure such as produced by function Surv of package survival.
#' @param v A (discretized) vector indexed as SS.
#' @param colors Colors to interpolate; must be a valid argument to col2rgb().
#' @param main A character string to explicit the title of the plot
#' @param legend A vector of character to explicit the legend of the plot
#' @param ... Parameters passed to plot function.
#' @importFrom survival survfit
#' @export
scurve = function(SS, v, colors=c("deepskyblue", "black", "red"), main="Survival", legend, ...) {
  pv = coxres(SS,v)[1]
  if (pv<1e-100) {
    pvt = "<1e-100"
  } else {
    pvt = format(pv,
    digits=3,scientific=TRUE)
  }
  sf=survfit(SS~v)
  levels = length(unique(v))
  col = colorRampPalette(colors)(levels)
  main= paste(main, " P=", pvt, sep="")
  plot(sf, col=col, main=main, ...)
  tab = table(v)
  if (missing(legend)) {
    if ("breaks" %in% names(attributes(v))) {
      b = signif(attr(v, "breaks"),3)
      legend = paste("[", b[1:(length(b)-1)], ",", b[2:length(b)], c(rep("[", length(b)-2), "]"), sep="")      
    } else {
      legend = names(tab)      
    }
  }
  legend = paste(legend, " (", tab, ")", sep="")
  legend("topright", legend=legend, col=col, pch=3, lty=1)
}

#' A Function That Fits the Cox Regression Model
#'
#' This function takes a survival structure such as produced by function Surv of
#' package survival. Takes a (discretized) vector v, indexed as SS.
#' Fits the Cox regression model of SS against v, and tests the
#' proportional hazards assumption.
#' Returns as a named vector of length 5:
#'      pvcox: the significance p-value (likelihood ratio test)
#'      pvhz: the p-value for the validity of the model
#'      hrlb: the lower bound of the 95% confidence interval for hr
#'      hr: the hazard ratio
#'      hrub: the upper bound of the 95% confidence interval for hr
#' @param SS A survival structure such as produced by function Surv of package survival.
#' @param v A (discretized) vector indexed as SS.
#' @return A named vector of length 5.
#' @importFrom survival coxph
#' @importFrom survival cox.zph
#' @export
coxres = function(SS,v) {
  f = suppressWarnings(coxph(SS~v))
  sf =  summary(f)
  pvcox = sf$logtest[3]
  tf = cox.zph(f)
  tf = tf$table
  pvhz = tf[dim(tf)[1],3]
  hr = sf$conf.int[1,1]
  hrlb = sf$conf.int[1,3]
  hrub = sf$conf.int[1,4]
  res = c(pvcox,pvhz,hrlb,hr,hrub)
  names(res) = c("pvcox","pvhz","hrlb","hr","hrub")
  return(res)
}

#' A Function Discretize a Vector of Numeric
#'
#' Takes a vector v of numeric. Replaces the values of v by integers
#' in (1:nd), according to values of breaks (increasing vector of numeric).
#' If breaks is missing, nd regular quantiles are taken.
#
#' @param v A vector of numeric.
#' @param nd The number of classes.
#' @param breaks An increasing vector of numeric.
#' @return A vector of factors.
#' @export
discr = function(v, nd=5, breaks){
  ind = !is.na(v)
  v1 = v[ind]
  if (missing(breaks)) {
    b = quantile(v1,(0:nd)/nd)
  } else {
    b = c(min(c(v1,breaks)), breaks, max(c(v1, breaks)))
  }
  b = unique(b)
  vd1 = cut(v1, b, include.lowest=TRUE, right=FALSE, labels=(1:(length(b)-1)))
  vd1 = as.vector(vd1)
  vd = v
  vd[ind] = vd1
  vd = as.factor(vd)
  return(structure(vd, breaks=b))
}

#' A Function That Returns the Longest Common Prefix.
#'
#' This function returns the longest common prefix.
#'
#' @param s A vector of character.
#' @param split character passed to strsplit function as split.
#' @param fixed boolean passed to gsub function as fixed.
#' @export
longest_common_prefix = function(s, split="", fixed=FALSE) {
  s = na.omit(s)
  l = min(sapply(strsplit(unique(s), split), length))
  f = lapply(strsplit(s, split), function(x) {x[1:l]})
  m = do.call(rbind, f)
  vl = apply(m, 2, function(col) {
    length(unique(col))
  })
  d_idx = which(vl!=1)
  if (length(d_idx) == 0) {
    d_idx = length(vl)
  }
  fd_idx = d_idx[1] - 1
  if (fd_idx > 0) {
    lcp = paste(c(f[[1]][1:fd_idx], ""), collapse=split)
  } else {
    lcp = ""
  }
  return(lcp)
}

#' A Function That Summarizes Expriment Grouping.
#'
#' This function summarizes expriment grouping.
#'
#' @param exp_grp a dataframe that describes the experimental grouping.
#' @export
summarize_exp_grp = function(exp_grp) {
  foo = sapply(colnames(exp_grp), function(cn) {
    nb_fact = length(unique(exp_grp[[cn]]))
    if (nb_fact > 1 & nb_fact < length(exp_grp[[cn]])) {
      print(paste("____________________________", cn, "____________________________"))
      print(table(exp_grp[[cn]]))
    }
  })
}
#' A Function That Simplifies Expriment Grouping.
#'
#' This function simplyfies expriment grouping.
#'
#' @param exp_grp a dataframe that describes the experimental grouping.
#' @return a simplified experimental grouping.
#' @export
simplify_exp_grp = function(exp_grp) {
  colnames(exp_grp) = simplify_column_names(colnames(exp_grp))
  col_div = apply(exp_grp, 2, function(col) {
    length(unique(col))
  })
  exp_grp = exp_grp[,col_div != 1]
  for (n in colnames(exp_grp)) {
    exp_grp[[n]] = simplify_column_names(exp_grp[[n]])

    exp_grp[[n]] = simplify_factor_names(exp_grp[[n]])
  }
  return(exp_grp)
}

#' A Function That Computes de False Discovery Rate.
#'
#' This function computes the false discovery rate from a vector of independent p-values.
#' It extracts the cutoff corresponding to the specified FDR. See Benjamini & Hochberg 1995 paper.
#' @param x A vector of independent p-values.
#' @param FDR The false discovery rate.
#' @return The corresponding cutoff.
#' @export
FDR = function (x, FDR)
{
    x <- sort(na.omit(x))
    N = length(x)
    i = 1
    while (N * x[i]/i < FDR & i <= N) i = i + 1
    if (i == 1)
        return(NA)
    else return(x[i - 1])
}

#' A Function That Builds a Fake Study.
#'
#' This function builds a fake study.
#'
#' @param nb_samples An integer that describes the number of samples.
#' @param nb_probes An integer that describes the number of probes.
#' @return a fake study.
#' @export
get_fake_study = function(nb_samples = 12, nb_probes = 10) {
  data = matrix(round(c(rnorm(nb_probes * floor(nb_samples/2)), rnorm(nb_probes * ceiling(nb_samples/2), 3,1)),3), nb_probes)
  data = data - min(data)
  colnames(data) = c(paste(rep("ctrl", nb_samples/2), 1:(nb_samples/2), sep=""), paste(rep("case", nb_samples/2), 1:(nb_samples/2), sep=""))
  rownames(data) = paste(rep("prb", nb_probes), 1:nb_probes, sep="")
  exp_grp = data.frame(
    sex=ifelse(runif(nb_samples)>0.5, "Male", "Female"),
    age=round(runif(nb_samples,20,30)),
    tabac=rep(c(rep("Smoker", nb_samples/4), rep("Non Smoker", nb_samples/4)), 2),
    treatment = c(rep("0 ug", nb_samples/2), rep("15 ug", nb_samples/2)),
    histo=rep("lung", nb_samples)
  )
  rownames(exp_grp) = colnames(data)
  platform = data.frame(
    gene_name = rep("...",nb_probes),
    GOID = rep("...",nb_probes)
  )
  rownames(platform) = rownames(data)
  return(list(data=data, exp_grp=exp_grp, platform=platform))
}
#' A Function That Computes `mean` + 2 * `sd` on a Numeric Vector.
#'
#' This function computes `mean` + 2 * `sd` on a numeric vector.
#'
#' @param ctrl A numeric vector
#' @return `mean` + 2 * `sd` of the inpuit vector
#' @export
m2sd = function(ctrl) {
  mean(ctrl) + 2 * sd(ctrl)
}
#' A Function That Simplifies Sample Names.
#'
#' This function simplyfies sample names.
#'
#' @param sample_names A character vector that describes the sample names.
#' @return A character vector that describes simplified sample names.
#' @export
simplify_sample_names = function(sample_names) {
  sample_names = sub(".CEL.gz", "", sample_names, ignore.case = TRUE)
  sample_names = sub(".CEL", "", sample_names, ignore.case = TRUE)
  tmp_sample_names = as.vector(sapply(sample_names, function(gsm) {
    as.list(strsplit(gsm, "_")[[1]][1])
  }))
  if (sum(duplicated(tmp_sample_names)) == 0) {
    sample_names = tmp_sample_names
  }
  # sample_names = do.call(cbind, t(sample_names))
  return(sample_names)
}
#' A Function That Simplifies Factor Names.
#'
#' This function simplyfies factor names.
#'
#' @param factor_names A character vector that describes the factor names.
#' @param split character passed to strsplit function as split.
#' @param fixed boolean passed to gsub function as fixed.
#' @return A character vector that describes simplified factor names.
#' @export
simplify_factor_names = function(factor_names, split="", fixed=FALSE) {
  factor_names = as.character(factor_names)
  one_way = function(factor_names, split, fixed) {
    # print(factor_names)
    l = min(sapply(strsplit(unique(factor_names), split), length))
    f = lapply(strsplit(factor_names, split), function(x) {x[1:l]})
    m = do.call(rbind, f)
    vl = apply(m, 2, function(col) {
      length(unique(col))
    })
    d_idx = which(vl!=1)
    if (length(d_idx) == 0) {
      d_idx = length(vl)
    }
    fd_idx = d_idx[1] - 1
    if (fd_idx > 0) {
      to_be_remove = paste(c(f[[1]][1:fd_idx], ""), collapse=split)
      if (!fixed) {
        to_be_remove = paste("^", to_be_remove, sep="")
      }
      factor_names = gsub(to_be_remove, "", factor_names, fixed=fixed)
      # factor_names = gsub(to_be_remove, "", factor_names, fixed=fixed)
    }
    # factor_names = substr(factor_names, d_idx, 100000L)
    # print(factor_names)
    return(factor_names)
  }
  str_rev = function(factor_names) {
    sapply(lapply(strsplit(factor_names, ""), rev), function(s) { paste(s, collapse='')})
  }
  factor_names = one_way(factor_names, split, fixed)
  factor_names = str_rev(factor_names)
  factor_names = one_way(factor_names, split, fixed)
  factor_names = str_rev(factor_names)
  # print(factor_names)

  return(factor_names)
}
#' A Function That Simplifies Column Names.
#'
#' This function simplyfies column names.
#'
#' @param column_names A character vector that describes the column names.
#' @return A character vector that describes simplified column names.
#' @export
simplify_column_names = function(column_names) {
  n = tolower(column_names)
  n = gsub("[.|:| ]", "_", n)
  n = gsub("_+", "_", n)
  n = gsub("^_", "", n)
  n = gsub("_$", "", n)
  return(n)
}
#' Retrieve Directory's .CEL.gz Files.
#'
#' This function takes a directory as parameter and return a string vector of the .CEL.gz files
#' that it contains.
#'
#' @param dir A character string to consider.
#' @param full.names A logical value.
#' @param ERROR_IF_EMPTY A logical value that allow to not throw an error if dir do not contain any .CEL.gz file.
#' @return It returns a string vector of the .CEL.gz files.
#' @export
get_cel_filenames = function(dir, full.names=TRUE, ERROR_IF_EMPTY=TRUE) {
  raw_dir = paste(dir, "/raw", sep="")
  cel_filenames = list.files(raw_dir, pattern="*.CEL.gz", full.names=full.names, ignore.case = TRUE)
  if (length(cel_filenames) == 0) {
    cel_filenames = list.files(dir, pattern="*.CEL.gz", full.names=full.names, ignore.case = TRUE)
  }
  if (ERROR_IF_EMPTY & length(cel_filenames) == 0 ) { stop(paste("ERROR! No .CEL.gz file in ", paste(dir, collapse=" "), ".")) }
  return(cel_filenames)
}

#' Fuse two experimental grouping
#'
#' This function tales two experimental grouping as parameter and merge them according to their row names.
#' @param exp_grp1 The first experimental grouping.
#' @param exp_grp2 The second experimental grouping.
#' @param by Specifications of the columns used for merging. See 'Details' in ?merge.
#' @return The fused experimental grouping
#' @export
fuse_exp_grp = function(exp_grp1, exp_grp2, by="row.names") {
  simplify_exp_grp_rownames = function(exp_grp2, exp_grp1) {
    rownames(exp_grp1) = sub(".CEL.gz", "", rownames(exp_grp1), ignore.case = TRUE)
    rownames(exp_grp2) = sub(".CEL.gz", "", rownames(exp_grp2), ignore.case = TRUE)
    new_names = sapply(rownames(exp_grp2), function(n) {
      rownames_grep = unlist(sapply(rownames(exp_grp1), grep, n))
      if (length(rownames_grep) == 1) {
        return(names(rownames_grep))
      } else if (length(rownames_grep) > 1) {
        stop(paste("More than one rowname match ", n, ".", sep=""))
      } else {
        return(n)
      }
    })
    return(new_names)
  }
  rownames(exp_grp1) = simplify_exp_grp_rownames(exp_grp1, exp_grp2)
  rownames(exp_grp2) = simplify_exp_grp_rownames(exp_grp2, exp_grp1)

  fused_exp_grp = merge(exp_grp1, exp_grp2, by=by, all=TRUE)
  rownames(fused_exp_grp) = fused_exp_grp[,"Row.names"]
  fused_exp_grp[,"Row.names"] = NULL
  return(fused_exp_grp)
}

#' Retrieving .CEL.gz Files from GEO.
#'
#' This function retrieves .CEL.gz files from GEO to the file system. It takes as parameters a character vector of GMSs to be retrived and a targeted directory. It returns a character vector of .CEL.gz files.
#' @param gsms A character vector of GMSs to be retrived.
#' @param basedir A string describing the targeted directory.
#' @return A character vector of .CEL.gz files.
#' @export
#' @importFrom GEOquery getGEOSuppFiles
get_gsm = function(gsms, basedir) {
  "Retrieve .CEL.gz from NCBI GEO web site using GSM identifiers."
  dir.create(path=basedir, showWarnings=FALSE, recursive=TRUE)
  cel_filenames = get_cel_filenames(basedir, ERROR_IF_EMPTY=FALSE)
  sapply(gsms, function(gsm)  {
    gsm_grep = grep(gsm, cel_filenames)
    if (length(gsm_grep) == 0) {
      getGEOSuppFiles(gsm, baseDir=basedir, makeDirectory=FALSE)
    } else if (length(gsm_grep) != 1) {
      stop(paste("More than one .CEL.gz file match ", gsm, " in ", basedir, sep=""))
    } else {
      print(paste("Using locally cached version: ", cel_filenames[gsm_grep], " for ", gsm, ".", sep=""))
    }
  })
  return(get_cel_filenames(basedir))
}

#' Broken-stick Distribution for 'p' Pieces.
#'
#' Compute the expected values of the broken-stick distribution for `p` pieces.
#
#  @author Pierre Legendre
#' @param p An interger that specify the number of pieces
#' @return broken-stick distribution for `p` pieces.
#' @examples
#' broken_stick_out_20 = broken_stick(20)
#' @export
broken_stick <- function(p)
{
result = matrix(0,p,2)
colnames(result) = c("j","E(j)")
for(j in 1:p) {
   E = 0
   for(x in j:p) E = E+(1/x)
   result[j,1] = j
   result[j,2] = E/p
   }
return(result)
}

#' A memoised version of getGEO
#'
#' This function offer a memoised version of getGEO. It is used to get GEO data
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to getGEO function.
#' @return GEO content.
#' @examples
#' # foo = mgetGEO("GSE42707", getGPL=FALSE)
#' @importFrom GEOquery getGEO
#' @importFrom memoise memoise
mgetGEO = memoise(function(...) {
  getGEO(...)
})

#' An other memoised version of getGEO
#'
#' This function offer a memoised version of getGEO. It returns a Table
#' version of a GEO query response. It is used to get platform metadatas
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to getGEO function.
#' @return GEO content.
#' @examples
#' # foo = mgetGEO("GPL13534", getGPL=FALSE)
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @importFrom memoise memoise
mtgetGEO = memoise(function(...) {
  Table(getGEO(...))
})


#' A Monitored Version of `apply`
#'
#' This function offer a monitored version of `apply`.
#'
#' @param mat A matrix usualy gives as `apply` first parameter.
#' @param marg An integer describing the marginal to use, usualy gives as  `apply` second parameter.
#' @param func A function usualy gives as `sapply` third parameter.
#' @param mod An integer that define the frequency of the monitoring.
#' @param ... Parameters passed to `func` function.
#' @return A vector of the application of the `func` function to each element of the `vec` vector.
# ' @examples
# ' # foo = monitored_apply(matrix(rnorm(50), 10), 1, function(v) {print(length(v)); Sys.sleep(1); return(c(mean(v), sd(v)))}, mod = 1)
#' @export
monitored_apply = function(mat, marg=1, func, mod=100, ...) {
  if (marg==1) {
    nb_it = nrow(mat)
    tmp_matrix = cbind(1:nb_it, mat)
  } else {
    nb_it = ncol(mat)
    tmp_matrix = rbind(1:nb_it, mat)
  }
  d1 = as.numeric(Sys.time())
  foo = apply(tmp_matrix, marg, function(vec) {
    id = as.numeric(vec[1])
    ret = func(vec[-1], ...)
    # monitoring...
    if (id %% mod == 0) {
      d2 = Sys.time()
      d2_num = as.numeric(Sys.time())
      elapse = d2_num - d1
      remain = (nb_it - id) * elapse / id
      end_at = d2 + remain
      print(paste("~", round(remain), " seconds remaining, finishing at ~", end_at, " (", id, "/", nb_it, ")." , sep=""))
    }
    return(ret)
  })
}

#' A Monitored Version of `sapply``
#'
#' This function offer a monitored version of `sapply`.
#'
#' @param vec A vector usualy gives as `sapply` first parameter.
#' @param func A function usualy gives as `sapply` second parameter.
#' @param ... Parameters passed to `monitored_apply` function.
#' @return A vector of the application of the `func` function to each element of the `vec` vector.
# ' @examples
# ' # foo = monitored_sapply(1:10, function(i) {print(i); Sys.sleep(1); return(c(i*i, i+i))}, mod = 1)
#' @export
monitored_sapply = function(vec, func, ...) {
  return(monitored_apply(t(vec), 2, func, ...))
}

#' A memoised version of readRDS
#'
#' This function offer a memoised version of readRDS. It is used to get GEO data
#' from GEO portal or a local file.
#'
#' @param ... parameters passed to readRDS function.
#' @return RDS content.
#' @examples
#' # foo = mreadRDS("GSE42707", getGPL=FALSE)
#' @importFrom memoise memoise
mreadRDS = memoise(function(...) {
  readRDS(...)
})