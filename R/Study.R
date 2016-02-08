
#' A memoised version of getGEO
#'

#' This fucntion offer a memoised version of getGEO. It is used to get GEO data 
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
#' This fucntion offer a memoised version of getGEO. It returns a Table 
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

#' Study Factory
#'
#' Intanciates an object that extends the epimedtools internal class Study_abstract.
#'
#' This function intanciates an object that extends the epimedtools internal abstract class Study_abstract. You can learn more about Study_abstract class by typing ?epimedtools::Study_abstract
#' 
#' @param cache_filename A character string that describes the study cache file name.
#' @param Study_RC_name A character string that specify the RC to use to instanciate the object.
#' @examples
#' # create study and load data from GSE57831_series_matrix.txt.gz contained in the package
#' study = create_study()
#' study$series_matrix_filename = system.file(
#'    "extdata/GSE26471", 
#'    "GSE26471_series_matrix.txt.gz", 
#'    package = "epimedtools"
#'  )  
#' # str(study$get_gset(dest_dir="data"))
#' 
#' # load data from GEO
#' study = create_study()
#' study$gse = "GSE26471" # 1 sample, expression
#' # study$gse = "GSE42707" # 1 sample
#' # study$gse = "GSE56382" # 2 samples
#' # study$gse = "GSE49585" # 2 samples
#' # str(study$get_gset(dest_dir="data"))
#' 
#' # cache study
#' study$cache_it("/tmp/tmp_cached_study.rds")
#' # load a cached study
#' study = create_study("/tmp/tmp_cached_study.rds")
#' @return It returns an object that extends the Study_abstract class.
#' @export
create_study = function(cache_filename, Study_RC_name="Study_geo") {
  if (missing(cache_filename)) {
    study = get(Study_RC_name)()
  } else {
    study = get(Study_RC_name)(cache_filename)
  }
  return(study)
}


#' A Reference Class to represent a multi-omic study.
#'
#' This class allow to manipulate useful concepts involved in a multi-omics
#' study. It an abstract class.
#' 
#' @field plaform The Table version of the GEO platform description.
#' @field platform_filenames A character string that describes the platform file name.
#' @field cache_filename A character string that describes the study cache file name.
#' @export
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @import methods
Study_abstract = setRefClass(
  "Study_abstract",
  fields = list(
    "platform" = "ANY",
    "platform_filename" = "character",
    "cache_filename" = "character"
  ),
  methods = list(
    initialize = function(cfn) {
      "Constructor."
      if (!missing(cfn)) {
        .self$cache_filename = cfn
        if (!file.exists(.self$cache_filename)) {
          cache_dir = paste(rev(rev(
            strsplit(.self$cache_filename, "/")[[1]]
          )[-1]), collapse = "/")
          dir.create(cache_dir, showWarnings = FALSE, recursive =
                       TRUE)
          cache_it()
        } else {
          print("reifiyng study...")
          s2 = readRDS(.self$cache_filename)
          # print(.self)
          for (f in names(s2$getRefClass()$fields())) {
            .self[[f]] = s2[[f]]
          }
          .self$cache_filename = cfn
          print("done.")
        }
      }
    },
    cache_it = function(cache_filename) {
      "Writes the study on the disk."
      if (!missing(cache_filename)) {
        .self$cache_filename = cache_filename
      }
      if (length(.self$cache_filename) == 1) {
        print("caching study...")
        saveRDS(.self, .self$cache_filename)
        print("done.")
      }
    },
    get_platform = function(CACHE = FALSE, MEMOISE = TRUE, dest_dir = "data/platforms") {
      "Computes if not and returns the platform field."
      if (is.null(dim(.self$platform))) {
        if (length(.self$platform_filename) == 1) {
          if (file.exists(.self$platform_filename)) {
            print(
              paste(
                "Launching platform informations from ", .self$platform_filename, "...", sep =
                  ""
              )
            )
            if (MEMOISE) {
              .self$platform = mtgetGEO(filename = .self$platform_filename, getGPL = FALSE)
            } else {
              .self$platform = Table(getGEO(
                filename = .self$platform_filename, getGPL = FALSE
              ))
            }
            print("done.")
          } else {
            stop(paste("No", .self$platform_filename, "file."))
          }
        } else {
          print(paste(
            "Launching ", .self$get_platform_name(), " from GEO...", sep =
              ""
          ))
          dir.create(dest_dir, showWarnings = FALSE, recursive =
                       FALSE)
          if (MEMOISE) {
            .self$platform = mtgetGEO(.self$get_platform_name(), getGPL = FALSE, destdir =
                                        dest_dir)
          } else {
            .self$platform = Table(getGEO(
              .self$get_platform_name(), getGPL = FALSE, destdir = dest_dir
            ))
            .self$platform_filename = paste(dest_dir, "/", .self$get_platform_name(), ".soft", sep="")
          }
        }
        if (CACHE) {
          cache_it()
        }
      }
      return(.self$platform)
    },
    get_platform_name = function(...) {
      "Abstract method that gives access to the platform name."
      stop("Call to the abstract method get_platform_name of the class Study_abstract.")
    },
    get_exp_grp = function(...) {
      "Abstract method that gives access to the experimental grouping."
      stop("Call to the abstract method get_exp_grp of the class Study_abstract.")
    },
    get_data = function(...) {
      "Abstract method that gives access to the data."
      stop("Call to the abstract method get_data of the class Study_abstract.")
    },
    # Analysis method dealing with Study_abstract class intances...
    plot_qc = function(method = "boxplot", ...) {
      "Plot the quality control of the study."
      get(method)(t(na.omit(.self$get_data())) ~ colnames(.self$get_data()), las =
                    2, ...)
    },
    get_ratio = function(ctrl_sample_names, method = "mean", ...) {
      "Normalize data accross probes using method *method* according to a reference poulation describe by  *ctrl_sample_names* its sample names vector."
      if (missing(ctrl_sample_names)) {
        d = .self$get_data(...)
      } else {
        d = .self$get_data(...)[, ctrl_sample_names]
      }
      ctrl_op = apply(d, 1, get(method))        
      ratio = .self$get_data(...) / ctrl_op
      return(ratio)
    }
  )
)




#' A Reference Class to represent a multi-omic study.
#'
#' This class extends Study_abstract by encapsulation GEO dataset.
#' 
#' @field gset The GEO set of data as it is return by getGEO.
#' @field gse A character string that describes the GEO accession number.
#' @field series_matrix_filename A character string that describes the GEO serie matrix family file name.
#' @export
#' @importFrom GEOquery getGEO
#' @importFrom Biobase exprs
#' @import methods
Study_geo = setRefClass(
  "Study_geo",
  fields = list(
    "gset" = "ANY",
    "gse" = "character",
    "series_matrix_filename" = "character"
  ),
  contains = "Study_abstract",
  methods = list(
    get_gset = function(CACHE=TRUE, MEMOISE=FALSE, dest_dir="data") {
      "Computes (if not yet done) and returns the gset field."
      if (is.null(dim(.self$gset))) {
        if (length(.self$series_matrix_filename) == 1) {
          if (file.exists(.self$series_matrix_filename)) {
            print(paste(
              "Launching data from ", .self$series_matrix_filename, "...", sep = ""
            ))
            if (MEMOISE) {
              .self$gset = mgetGEO(filename = .self$series_matrix_filename, getGPL = FALSE)
            } else {
              .self$gset = getGEO(filename = .self$series_matrix_filename, getGPL = FALSE)
            }
            print("done.")
          } else {
            stop(paste("No", .self$series_matrix_filename, "file."))
          }
        } else if (length(.self$gse) == 1) {
          print(paste("Launching ", .self$gse, " from GEO...", sep = ""))
          dest_dir_gse = paste(dest_dir, "/", .self$gse, sep="")
          dir.create(dest_dir_gse, showWarnings = FALSE, recursive =
                       TRUE)
          if (MEMOISE) {
            tmp_gset = mgetGEO(.self$gse, getGPL = FALSE, destdir=dest_dir_gse)
          } else {
            tmp_gset = getGEO(.self$gse, getGPL = FALSE, destdir=dest_dir_gse)
          }
          if (length(.self$gset) == 1) {
            .self$gset = tmp_gset[[1]]
            .self$series_matrix_filename = paste(dest_dir_gse, "/", names(tmp_gset)[1], sep =
                                                   "")
          } else {
            stop(paste(.self$gse, " is a SuperSeries... Not yet supported.", sep = ""))
          }
          print(paste(
            "done. File locally cached here: ", .self$series_matrix_filename, sep =
              ""
          ))
        } else {
          stop("You need to define gse (length 1) field to get it from GEO")
        }
        if (CACHE) {
          .self$cache_it()
        }
      }
      return(.self$gset)
    },
    get_platform_name = function(...) {
      "Proxy to annotation hold in gset."
      return(.self$get_gset(...)@annotation)
    },
    get_exp_grp = function(...) {
      "Proxy to phenoData@data hold in gset."
      return(.self$get_gset(...)@phenoData@data)
    },
    get_data = function(...) {
      "Proxy to data hold in gset."
      return(exprs(.self$get_gset(...)))
    }
  )
)



#' A Reference Class to represent a multi-omic study.
#'
#' This class extends Study_abstract by encapsulation hand made fields: "data", "platform_name", "exp_grp".
#' 
#' @field data A matris of experiemental data. Each column is a sample, named by it's sample name, and each row is a probes.
#' @field platform_name A character string that describes the platform used to perform.
#' @field exp_grp A data.frame that describe the experimental grouping. Each line is a sample.
#' @export
Study_loc = setRefClass("Study_loc",
  fields = list(
    "data" = "ANY",
    "platform_name" = "character",
    "exp_grp" = "ANY"
  ),
  contains = "Study_abstract",
  methods = list( 
    get_data = function() {
      return(data)
    },
    get_platform_name = function() {
      return(platform_name)
    },
    get_exp_grp = function() {
      return(exp_grp)
    }
  )
)

