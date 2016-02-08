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

