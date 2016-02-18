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
#' study$save("/tmp/tmp_cached_study.rds")
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
          .self$save()
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
    save = function(cache_filename) {
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
          .self$save()
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
    },
    plot_boxplot = function(probe_name, factor_name, ylim, bp_function_name="boxplot") {
      "Draw the box plot for a given `probe_name` and a given `factor_name`."
    	idx_sample = rownames(.self$get_exp_grp())
    	# Dealing with ratio data, reduce data to interesting probes and samples
    	filtred_bp_data = .self$get_data()[probe_name, idx_sample]
    	# Box plots
      if (missing(ylim)) {
    	  ylim = c(min(filtred_bp_data), max(filtred_bp_data))
      }
	    gene = "genename"
	    # boxplot_filename <- paste(study_dirname, "/", gene, "_", probe, "_", factor_name, ".pdf", sep="")
	    # pdf(file=boxplot_filename, height=10, width=10)# open jpeg device with specified dimensions
      bp_function = get(bp_function_name)
	    bp_function(unlist(filtred_bp_data[probe,])~get_exp_grp()[,factor_name], # open boxplot
          # what=c(1,1,1,0),
	      ylim = ylim,                  # fixed vertical scale
	      col = "grey", border = "black",  # colors of inside and border of box
	      las = 2,                    # written vertically
	      xlab = factor_name,
	      main = paste(gene, probe, sep="@") # title
	    )
	    # dev.off()
    }  
  )
)




