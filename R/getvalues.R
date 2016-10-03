#' @include getvalues.R
NULL

#' Evaluate an \code{rseHandle} object and return a
#' \code{RangedSummarizedExperiment}
#'
#' \code{getvalues} evaluates any valid  \code{rseHandle} object
#' and loads data into memory that is ordinarily located on disk.
#' 
#' @param handle An \code{rseHandle} object. 
#' @param format = "sparse" or "normal". Representation of data
#' when evaluted. 
#'
#' @return Returns a \code{RangedSummarizedExperiment} object. 
#'
#' @export
setGeneric(name = "getvalues", def = function(handle, format = "sparse")
    standardGeneric("getvalues"))

#' @rdname getvalues
setMethod("getvalues", c("rseHandle", "character"),
        definition = function(handle, format = "sparse") {
            return(NULL)
            #con <- dbConnect(drv=RSQLite::SQLite(), dbname=dblu)
            #sqlcmd <- paste0('select * from ', tblName,' where row in (', paste(rows, collapse=','), ') and column in (', paste(columns, collapse=','),')')
            #a <- dbGetQuery(conn=con, statement=sqlcmd)
            #dbDisconnect(con)
})