#' @include constructor.R
NULL

#' Subset an \code{rseHandle} object
#'
#' \code{rseHandleSubset} subsets the \code{rowData} or \code{colData}
#' based on the indicies of i (rowData) and j (colData). Can
#' also use the `[` operator as is common in R object subsetting.
#' 
#' @param x An \code{rseHandle} object to be subset. 
#' @param i Numeric indices to be subset from rows. 
#' @param j Numeric indices to be subset from columns
#' @param drop Other argument to pass to [
#'
#' @return Returns a subsetted \code{sqlSparseHandle} object. 
#' 
#' @examples
#' dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
#' d <- dat2[1:5, 3:4]
#' d <- rseHandleSubset(dat2, 1:5, 3:4)
#'
#' @export
setGeneric(name = "rseHandleSubset", def = function(x, i, j)
    standardGeneric("rseHandleSubset"))

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "numeric", "numeric"),
        definition = function(x, i, j) {

        existing <- as(x, "RangedSummarizedExperiment")
        if(i[1] == -999) i <- 1:dim(x)[1]
        if(j[1] == -999) j <- 1:dim(x)[2]
        rse <- existing[i,j]
        
        # Update mappings
        new_rowMap <- x@rowMap[as.character(i)]
        names(new_rowMap) <- as.character(1:length(i))
        new_colMap <- x@colMap[as.character(j)]
        names(new_colMap) <- as.character(1:length(j))
        
        b <- new("rseHandle", rse, rowMap = new_rowMap, colMap = new_colMap)
        return(b)
})

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "missing", "numeric"),
        definition = function(x, i, j) {
        return(rseHandleSubset(x,i=-999,j))
})

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "numeric", "missing"),
        definition = function(x, i, j) {
        return(rseHandleSubset(x,i,j=-999))
})

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "missing", "logical"),
        definition = function(x, i, j) {
        return(rseHandleSubset(x,i=-999,which(j)))
})

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "logical", "missing"),
        definition = function(x, i, j) {
        return(rseHandleSubset(x,which(i),j=-999))
})

#' @rdname rseHandleSubset
setMethod("rseHandleSubset", signature("rseHandle", "logical", "logical"),
        definition = function(x, i, j) {
        return(rseHandleSubset(x,which(i),which(j)))
})

#' @rdname rseHandleSubset
setMethod("[", signature(x = "rseHandle", i = "numeric", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,i,j))
})

#' @rdname rseHandleSubset
setMethod("[", signature(x = "rseHandle", i = "missing", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,i=-999,j))
})

#' @rdname rseHandleSubset 
setMethod("[", signature(x = "rseHandle", i = "numeric", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,i,j=-999))
})

#' @rdname rseHandleSubset
setMethod("[", signature(x = "rseHandle", i = "logical", j = "logical", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,which(i),which(j)))
})

#' @rdname rseHandleSubset
setMethod("[", signature(x = "rseHandle", i = "missing", j = "logical", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,i=-999,which(j)))
})

#' @rdname rseHandleSubset 
setMethod("[", signature(x = "rseHandle", i = "logical", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(rseHandleSubset(x,which(i),j=-999))
})


######################
# Subsetting by overlaps
######################

#' Subset an \code{rseHandle} object by region defined in a 
#' \code{GRanges} object
#' 
#' @param query An \code{rseHandle} object
#' @param subject A \code{GRanges} object of the regions of intrest
#' @param maxgap Refer to \code{GRanges} documentation of this function. 
#' @param minoverlap Refer to \code{GRanges} documentation of this function. 
#' @param type Refer to \code{GRanges} documentation of this function. 
#' @param ... Refer to \code{GRanges} documentation of this function. 
#' 
#' @return A \code{rseHandle} object subsetted
#' after the subjet was filtered 
#' @importMethodsFrom IRanges findOverlaps
#' @importMethodsFrom IRanges subsetByOverlaps
#' 
setMethod("subsetByOverlaps", signature = c("rseHandle", "GRanges"),
    function(query, subject, ...) {
        
        existing <- as(query, "RangedSummarizedExperiment")
        rse <- subsetByOverlaps(existing, subject, ...)
        mmap <- findOverlaps(existing, subject, ...)
            
        # Update mappings
        i <- from(mmap)
        new_rowMap <- query@rowMap[as.character(i)]
        names(new_rowMap) <- as.character(1:length(i))
        
        b <- new("rseHandle", rse, rowMap = new_rowMap, colMap = query@colMap)
        return(b)
})
