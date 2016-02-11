#' Return directory's .CEL.gz files.
#'
#' This method takes a directory as parameter and return a string vector of the .CEL.gz files 
#' that it contains.
#'
#' @param dir A character string to consider.
#' @param full.names A logical value. 
#' @return It returns a string vector of the .CEL.gz files. 
#' @export
get_cel_filenames = function(dir, full.names=TRUE) {
  raw_dir = paste(dir, "/raw", sep="")
  cel_filenames = list.files(raw_dir, pattern="*.CEL.gz", full.names=full.names)
  if (length(cel_filenames) == 0) {
    cel_filenames = list.files(dir, pattern="*.CEL.gz", full.names=full.names)
  }
  if (length(cel_filenames) == 0 ) { stop(paste("ERROR! No .CEL.gz file in ", paste(dir, collapse=" "), ".")) }
  return(cel_filenames)
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
    set_data = function(celfile.path=getwd()) {
      "Set the data field and update experiment grouping."
      cel_files = lapply(.self$cel_filedirs, get_cel_filenames)
      orig = unlist(sapply(1:length(.self$cel_filedirs), function(i) {
        split1 = unique(unlist(strsplit(.self$cel_filedirs[i], "/")))  
        split2 = unique(unlist(strsplit(.self$cel_filedirs[-i], "/")))  
        n = paste(split1[!(split1 %in% split2)], collapse="_")
        rep(n, length(cel_files[[i]]))
      }))
      cel_files = unlist(cel_files)
      cel_files_short = unlist(lapply(.self$cel_filedirs, get_cel_filenames, full.names=FALSE))
      cel_files = cel_files[!duplicated(cel_files_short)]
      orig = orig[!duplicated(cel_files_short)]
      sample_names = cel_files_short[!duplicated(cel_files_short)]
      tmp_exp_grp = data.frame(orig=orig)
      rownames(tmp_exp_grp) = sample_names 
      .self$exp_grp = tmp_exp_grp
      .self$data = exprs(justRMA(filenames=cel_files, celfile.path=celfile.path))
    }
  )
)


