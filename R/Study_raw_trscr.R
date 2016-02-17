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
  cel_filenames = list.files(raw_dir, pattern="*.CEL.gz", full.names=full.names)
  if (length(cel_filenames) == 0) {
    cel_filenames = list.files(dir, pattern="*.CEL.gz", full.names=full.names)
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
  fused_exp_grp = merge(exp_grp1, exp_grp2, by=by, all=TRUE)
  rownames(fused_exp_grp) = fused_exp_grp[,"Row.names"]
  fused_exp_grp[,"Row.names"] = NULL  
  return(fused_exp_grp)
}


#' A Reference Class to represent a multi-omic study.
#'
#' This class extends Study_loc by offering methods allowing to build data from cel files.
#' 
#' @field cel_filedirs A character string that describes the directory taht contain associated cel files.
#' @export
#' @importFrom Biobase exprs
#' @importFrom affy justRMA
Study_raw_trscr = setRefClass("Study_raw_trscr",
  fields = list(
    cel_filedirs="character"
  ),
  contains = "Study_loc",
  methods = list( 
    set_data = function() {
      "Set the data field and update experiment grouping."
      cel_files = lapply(.self$cel_filedirs, get_cel_filenames)
      orig = unlist(sapply(1:length(.self$cel_filedirs), function(i) {
        split1 = unique(unlist(strsplit(.self$cel_filedirs[i], "/")))  
        split2 = unique(unlist(strsplit(.self$cel_filedirs[-i], "/")))  
        n = paste(split1[!(split1 %in% split2)], collapse="_")
        rep(n, length(cel_files[[i]]))
      }))
      cel_files = unlist(cel_files)
      cel_files = sapply(cel_files, function(cel_file) {
        if (substr(cel_file, 1, 1) != "/") {
          return(paste(getwd(), "/", cel_file, sep=""))
        } else {
          return(cel_file)
        }
      })
      cel_files_short = unlist(lapply(.self$cel_filedirs, get_cel_filenames, full.names=FALSE))
      cel_files = cel_files[!duplicated(cel_files_short)]
      orig = orig[!duplicated(cel_files_short)]
      sample_names = cel_files_short[!duplicated(cel_files_short)]
      tmp_exp_grp = data.frame(orig=orig)
      rownames(tmp_exp_grp) = sample_names 
      if (!is.null(dim(.self$get_exp_grp()))) {
        tmp_exp_grp = fuse_exp_grp(.self$get_exp_grp(), tmp_exp_grp)
      }
      .self$exp_grp = tmp_exp_grp
      .self$data = exprs(justRMA(filenames=cel_files, celfile.path=""))
    }
  )
)


