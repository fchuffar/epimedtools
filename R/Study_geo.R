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
          .self$save()
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

