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
  data = matrix(round(rnorm(nb_probes * nb_samples),3), nb_probes)
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
#' A Function That Simplify Sample Names.
#'
#' This function simplfy sample names.
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
#' This function retrieves .CEL.gz files from GEO to the file system. It takes as parameters a character vector of GMSs to be retrived and a targeted directory. It retruns a character vector of .CEL.gz files.
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