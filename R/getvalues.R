#' @include getvalues.R
NULL

# Grab column and rows from file specified in rn character vector
# 3 is the type; 4 is the format; 1 is the file name; 2 is the index/table
# Sample names are needed for sqlite normal matrix
getWrapper <- function(row, column, rn, sampleNames){
    if(rn[3] == "sparse" & rn[4] == "sqlite"){
        con <- dbConnect(drv=RSQLite::SQLite(), dbname=rn[1])
        sqlcmd <- paste0('select * from ', rn[2],' where row in (', paste(row, collapse=','), ') and column in (', paste(column, collapse=','),')')
        mat <- dbGetQuery(conn=con, statement=sqlcmd)
        dbDisconnect(con)
        return(mat)
    } else if(rn[3] == "normal" & rn[4] == "HDF5"){
        return(h5read(rn[1], rn[2], index=list(row,column)))
    } else if(rn[3] == "normal" & rn[4] == "sqlite"){
        samples <- sampleNames[column]
        con <- dbConnect(drv=RSQLite::SQLite(), dbname=rn[1])
        sqlcmd <- paste0('select ', paste(samples, collapse=','),' from ', rn[2],' where row_names in (', paste(row, collapse=','), ')') 
        mat <- dbGetQuery(conn=con, statement=sqlcmd)
        dbDisconnect(con)
        return(mat)
    } else {
        stop("Invalid backend file specification")   
    }
}

cleanup <- function(dat, cols, rowMap, colMap, type, format){
    s_rowMap <- as.numeric(names(rowMap))
    names(s_rowMap) <- as.character(unname(rowMap))
    s_colMap <- as.numeric(names(colMap))
    names(s_colMap) <- as.character(unname(colMap))
    
    sparse_matrices <- lapply(1:length(type), function(i){
        samples <- cols[[i]]
        if(type[i] == "sparse"){
            scm <- s_colMap[samples]
            lmat <- cbind(s_rowMap[as.character(dat[[i]]$row)], scm[as.character(dat[[i]]$column)], dat[[i]]$value)
        } else { # normal
            d <- dat[[i]]
            colnames(d) <- samples
            rownames(d) <- seq(1, dim(d)[1], 1)
            md <- reshape2::melt(data.matrix(d))
            lmat <- md[md[,3] != 0,]
            colnames(lmat) <- NULL
        }
            rownames(lmat) <- NULL
            return(data.matrix(lmat))
    })
    sm <- do.call("rbind", sparse_matrices)
        
    if(format == "sparse"){
        return(Matrix::sparseMatrix(i = sm[,1], j = sm[,2], x = sm[,3]))
    } else {
        return(data.matrix(Matrix::sparseMatrix(i = sm[,1], j = sm[,2], x = sm[,3])))
    }
}


#' Evaluate an \code{rseHandle} object and return a
#' \code{RangedSummarizedExperiment} or a \code{seHandle}
#' object and rturn a \code{SummarizedExperiment}
#'
#' \code{getvalues} evaluates any valid  \code{yoloHandle} object
#' and loads data into memory that is ordinarily located on disk.
#' 
#' @param handle An  \code{seHandle} or \code{rseHandle} object. 
#' @param format = "sparse" or "normal". Representation of data
#' when evaluted. 
#'
#' @return Returns a \code{(Ranged)SummarizedExperiment} object. 
#'
#' @importFrom reshape2 acast melt
#' @importFrom rhdf5 h5read
#'
#' @export
setGeneric(name = "getvalues", def = function(handle, format = "sparse")
    standardGeneric("getvalues"))

#' @rdname getvalues
setMethod("getvalues", c("yoloHandle", "ANY"),
        definition = function(handle, format = "sparse") {
            
            beast <- ifelse(as.character(class(handle)) == "rseHandle", "rse", "se")
            stopifnot(format %in% c("sparse", "normal"))
            holyfour <- c("lookupFileName", "lookupTableName", "lookupFileType", "lookupFileFormat")
            
            # Define lookup table (lt) and iterate over row names (rn)
            lt <- data.frame(colData(handle)[,holyfour])
            rn <- unique(lt)
            rownames(rn) <- NULL
            n <- dim(rn)[1]
            
            # Get columns/lookup
            cols <- lapply(1:n, function(i){
                which(rn[i,1] == lt[,1] & rn[i,2] == lt[,2] & rn[i,3] == lt[,3] & rn[i,4] == lt[,4])
            })
            
            # Get data/lookup
            rowtrans <- unname(handle@rowMap)
            dat <- lapply(1:n, function(i){
                col <- cols[[i]]
                coltrans <- unname(handle@colMap[col])
                getWrapper(rowtrans, coltrans, as.character(rn[i,]), colnames(handle))
            })
            
            # Make final matrix with correct indicies; drop lookups; return RSE
            mat <- cleanup(dat = dat, cols = cols, rowMap = handle@rowMap, colMap = handle@colMap, type = rn$lookupFileType, format = format)
            if(beast == "rse"){
                x <- as(handle, "RangedSummarizedExperiment")
            } else {
                x <- as(handle, "SummarizedExperiment")
            }
            x@assays <- Assays(list(data=mat))
            x@colData <- x@colData[,-which(names(x@colData) %in% holyfour)]
            return(x)
})

