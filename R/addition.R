#' @include subsetting.R
NULL

#' Combining together multiple rseHandle objects
#' 
#' To add multiple \code{rseHandles} together, two rather stringent conditions
#' must be met. First, the \code{rowRanges} slots must be identical.
#' Moreover, the column names of the \code{colData} slot must also
#' be identical. The sames can then be merged into one object, 
#' the the varian
#' 
#' @param e1 A \code{rseHandles} object
#' @param e2 A \code{rseHandles} object
#' 
#' @return An rseHandle object joined
#' 
#' @examples 
#' dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
#' dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
#' dat4 <- dat1 + dat3
#' 
setMethod("+", signature(e1 = "rseHandle", e2 = "rseHandle"),
          definition = function(e1, e2) {
              
    #Basic Error Handling
    stopifnot(all(e1@rowRanges == e2@rowRanges))
    stopifnot(all(names(e1@colData) == names(e2@colData)))
              
    # Constructor
    cdn <- rbind(e1@colData, e2@colData)
    rse <- SummarizedExperiment(rowRanges=e1@rowRanges, colData=cdn)
    rowMap <- c(e1@rowMap)
    colMap <- c(e2@rowMap)
    b <- new("rseHandle", rse, rowMap = rowMap, colMap = colMap)
    return(b)
})