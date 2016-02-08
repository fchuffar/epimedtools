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

