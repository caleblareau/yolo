#' @include subsetting.R
NULL

#' Combining together multiple rseHandle objects
#' 
#' To add multiple \code{rseHandles} together, two rather stringent conditions
#' must be met. First, the \code{rowRanges} slots must be identical.
#' Moreover, the column names of the \code{colData} slot must also
#' be identical. The samples can then be merged into one object, 
#' where e2's samples will be added to the end. The only time that
#' order of operations for addition matters is that the meta data of
#' the first element being added (e1) will be retained if different
#' from the meta data in e2. 
#' 
#' @param e1 A \code{rseHandle} object
#' @param e2 A \code{rseHandle} object
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
    rse <- as(e1, "RangedSummarizedExperiment")
    rse@colData <- cdn 
    rowMap <- c(e1@rowMap) # should be the same either way
    colMap <- c(e1@colMap,e2@colMap)
    e2names <- as.character(as.numeric(names(e2@colMap))+length(e1@colMap))
    names(colMap) <- c(names(e1@colMap), e2names)
    b <- new("rseHandle", rse, rowMap = rowMap, colMap = colMap)
    return(b)
})
