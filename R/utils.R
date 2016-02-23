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
