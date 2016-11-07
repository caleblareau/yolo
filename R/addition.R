#' @include subsetting.R
NULL

#' Combining together multiple rseHandle objects
#' 
#' To add multiple \code{se/rseHandles} together, two rather stringent conditions
#' must be met. First, the \code{rowData/rowRanges} slots must be identical.
#' Moreover, the column names of the \code{colData} slot must also
#' be identical. The samples can then be merged into one object, 
#' where e2's samples will be added to the end. The only time that
#' order of operations for addition matters is that the meta data of
#' the first element being added (e1) will be retained if different
#' from the meta data in e2. 
#' 
#' @param e1 A \code{yoloHandle} object
#' @param e2 A \code{yoloHandle} object of the same type
#' 
#' @return An yoloHandle object joined by the column space
#' 
#' @examples 
#' dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
#' dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
#' dat4 <- dat1 + dat3
#' 
setMethod("+", signature(e1 = "yoloHandle", e2 = "yoloHandle"),
          definition = function(e1, e2) {
    
    beast1 <- ifelse(as.character(class(e1)) == "rseHandle", "rse", "se")
    beast2 <- ifelse(as.character(class(e2)) == "rseHandle", "rse", "se")
    if(beast1 != beast2) stop("Cannot add an rseHandle to an seHandle!")
    
    #Basic Error Handling
    if(beast1 == "rse") stopifnot(all(e1@rowRanges == e2@rowRanges))
    if(beast1 == "se") stopifnot(all(all(rowData(e1) == rowData(e2))))
    stopifnot(all(names(e1@colData) == names(e2@colData)))
              
    # Constructor
    if(beast1 == "rse"){
        cdn <- rbind(e1@colData, e2@colData)
        rse <- as(e1, "RangedSummarizedExperiment")
        rse@colData <- cdn 
        rowMap <- c(e1@rowMap) # should be the same either way
        colMap <- c(e1@colMap,e2@colMap)
        e2names <- as.character(as.numeric(names(e2@colMap))+length(e1@colMap))
        names(colMap) <- c(names(e1@colMap), e2names)
        b <- new("rseHandle", rse, rowMap = rowMap, colMap = colMap)
    } else {
        cdn <- rbind(e1@colData, e2@colData)
        se <- as(e1, "SummarizedExperiment")
        se@colData <- cdn
        rowMap <- c(e1@rowMap) # should be the same either way
        colMap <- c(e1@colMap,e2@colMap)
        e2names <- as.character(as.numeric(names(e2@colMap))+length(e1@colMap))
        names(colMap) <- c(names(e1@colMap), e2names)
        b <- new("seHandle", se, rowMap = rowMap, colMap = colMap)
    }
    
    return(b)
})
