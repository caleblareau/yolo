#' @include constructor.R
NULL

#' Subset a \code{yoloHandle} object
#'
#' \code{yoloHandleSubset} subsets the \code{rowData} or \code{colData}
#' based on the indicies of i (rowData) and j (colData). Can
#' also use the `[` operator as is common in R object subsetting.
#' 
#' @param x A \code{yoloHandle} object to be subset. 
#' @param i Numeric indices to be subset from rows. 
#' @param j Numeric indices to be subset from columns
#' @param drop Other argument to pass to [
#'
#' @return Returns a subsetted \code{yoloHandle} object. 
#' 
#' @examples
#' dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
#' d <- dat2[1:5, 3:4]
#' d <- yoloHandleSubset(dat2, 1:5, 3:4)
#'
#' @export
setGeneric(name = "yoloHandleSubset", def = function(x, i, j)
    standardGeneric("yoloHandleSubset"))

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "numeric", "numeric"),
        definition = function(x, i, j) {
        
        beast <- ifelse(as.character(class(x)) == "rseHandle", "rse", "se")
        if(beast == "rse"){
            existing <- as(x, "RangedSummarizedExperiment")
        } else {
            existing <- as(x, "SummarizedExperiment")
        }
        
        if(i[1] == -999) i <- 1:dim(x)[1]
        if(j[1] == -999) j <- 1:dim(x)[2]
        xx <- existing[i,j]
        
        # Update mappings
        new_rowMap <- x@rowMap[as.character(i)]
        names(new_rowMap) <- as.character(1:length(i))
        new_colMap <- x@colMap[as.character(j)]
        names(new_colMap) <- as.character(1:length(j))
        
        # Return subsetted value
        if(beast == "rse"){
            b <- new("rseHandle", xx, rowMap = new_rowMap, colMap = new_colMap)
        } else {
            b <- new("seHandle", xx, rowMap = new_rowMap, colMap = new_colMap)
        }
        return(b)
})

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "missing", "numeric"),
        definition = function(x, i, j) {
        return(yoloHandleSubset(x,i=-999,j))
})

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "numeric", "missing"),
        definition = function(x, i, j) {
        return(yoloHandleSubset(x,i,j=-999))
})

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "missing", "logical"),
        definition = function(x, i, j) {
        return(yoloHandleSubset(x,i=-999,which(j)))
})

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "logical", "missing"),
        definition = function(x, i, j) {
        return(yoloHandleSubset(x,which(i),j=-999))
})

#' @rdname yoloHandleSubset
setMethod("yoloHandleSubset", signature("yoloHandle", "logical", "logical"),
        definition = function(x, i, j) {
        return(yoloHandleSubset(x,which(i),which(j)))
})

#' @rdname yoloHandleSubset
setMethod("[", signature(x = "yoloHandle", i = "numeric", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,i,j))
})

#' @rdname yoloHandleSubset
setMethod("[", signature(x = "yoloHandle", i = "missing", j = "numeric", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,i=-999,j))
})

#' @rdname yoloHandleSubset 
setMethod("[", signature(x = "yoloHandle", i = "numeric", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,i,j=-999))
})

#' @rdname yoloHandleSubset
setMethod("[", signature(x = "yoloHandle", i = "logical", j = "logical", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,which(i),which(j)))
})

#' @rdname yoloHandleSubset
setMethod("[", signature(x = "yoloHandle", i = "missing", j = "logical", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,i=-999,which(j)))
})

#' @rdname yoloHandleSubset 
setMethod("[", signature(x = "yoloHandle", i = "logical", j = "missing", 
    drop = "missing"), definition = function(x, i, j, drop) {
        return(yoloHandleSubset(x,which(i),j=-999))
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
#' @param invert Refer to \code{GRanges} documentation of this function. 
#' @param ... Refer to \code{GRanges} documentation of this function. 
#' 
#' @return A \code{rseHandle} object subsetted
#' after the subjet was filtered 
#' @importMethodsFrom IRanges findOverlaps
#' @importMethodsFrom IRanges subsetByOverlaps
#' 
#' @examples 
#' library(GenomicRanges)
#' chr12reg <- GRanges(seqnames=c("chr12"),
#'      ranges = IRanges(start=c(56150000),end=c(56500000)))
#' dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
#' ss <- subsetByOverlaps(dat2, chr12reg)
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
