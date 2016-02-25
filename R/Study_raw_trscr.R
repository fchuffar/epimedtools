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
  contains = "Study_abstract",
  methods = list( 
    set_data = function(hook, ...) {
      "Set the data field and update experiment grouping."
      if (!missing(hook)) {
        get(hook)(...)
      }
      cel_files = lapply(.self$cel_filedirs, get_cel_filenames)
      orig = unlist(sapply(1:length(.self$cel_filedirs), function(i) {
        split1 = unique(unlist(strsplit(.self$cel_filedirs[i], "/")))  
        split2 = unique(unlist(strsplit(.self$cel_filedirs[-i], "/")))  
        n = paste(split1[!(split1 %in% split2)], collapse="_")
        n = gsub("[.]", "_", n)
        n = gsub("__*", "_", n)
        n = gsub("^_", "", n)
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
      rownames(tmp_exp_grp) = simplify_sample_names(sample_names) 
      if (!is.null(dim(.self$get_exp_grp()))) {
        tmp_exp_grp = fuse_exp_grp(.self$get_exp_grp(), tmp_exp_grp)
      }
      .self$exp_grp = tmp_exp_grp
      .self$data = exprs(justRMA(filenames=cel_files, celfile.path=""))
      colnames(.self$data) = simplify_sample_names(colnames(.self$data))
    }
  )
)
