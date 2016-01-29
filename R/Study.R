
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
mgetGEO = memoise(getGEO)

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
#' Intanciates an object of the epimedtools internal class Study.
#'
#' This function intanciates an object of the epimedtools internal class Study. You can learn more about Study class by typing ?epimedtools::Study
#' 
#' @param cache_filename A character string that describes the study cache file name.
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
#' @return It returns an object of the Study class.
#' @export
create_study = function(cache_filename) {
  if (missing(cache_filename)) {
    study = Study()
  } else {
    study = Study(cache_filename)
  }
  return(study)
}


#' A Reference Class to represent a multi-omic study.
#'
#' This class allow to manipulate useful concepts involved in a multi-omics
#' study.
#' 
#' @field gset The GEO set of data as it is return by getGEO.
#' @field plaform The Table version of the GEO platform description.
#' @field gse A character string that describes the GEO accession number.
#' @field platform_filenames A character string that describes the platform file name.
#' @field series_matrix_filename A character string that describes the GEO serie matrix family file name.
#' @field cache_filename A character string that describes the study cache file name.
#' @export
#' @importFrom GEOquery Table
#' @importFrom GEOquery getGEO
#' @importFrom Biobase exprs
#' @import methods
Study = setRefClass(
  "Study",
  fields = list(
    "gset" = "ANY",
    "platform" = "ANY",
    "gse" = "character",
    "platform_filename" = "character",
    "series_matrix_filename" = "character",
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
          for (f in names(.self$getRefClass()$fields())) {
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
    get_gset = function(CACHE=TRUE, MEMOISE=FALSE, dest_dir="../../data") {
      "Computes (if not yet done) and returns the gset field."
      if (is.null(dim(.self$gset))) {
        if (length(.self$series_matrix_filename) == 1) {
          if (file.exists(.self$series_matrix_filename)) {
            print(paste(
              "Launching data form ", .self$series_matrix_filename, "...", sep = ""
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
          print(paste("Launching ", .self$gse, " form GEO...", sep = ""))
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
    get_platform = function(CACHE = FALSE, MEMOISE = TRUE, dest_dir = "../../data/platforms") {
      "Computes if not and returns the platform field."
      if (is.null(dim(.self$platform))) {
        if (length(.self$platform_filename) == 1) {
          if (file.exists(.self$platform_filename)) {
            print(
              paste(
                "Launching platform informations form ", .self$platform_filename, "...", sep =
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
            "Launching ", .self$get_platform_name(CACHE = FALSE), " form GEO...", sep =
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
    },
    plot_qc = function(method = "boxplot", ...) {
      "Plot the quality control of the study."
      get(method)(t(na.omit(.self$get_data())) ~ colnames(.self$get_data()), las =
                    2, ...)
    }
  )
)
