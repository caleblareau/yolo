#' @include classDefinitions.R
NULL

#' Create a \code{sqlSparseHandle} object
#'
#' \code{sqlSparseHandleMake} takes the three parameters specified below
#' and creates a \code{sqlSparseHandle} object which can be subset and
#' manipulated without reading any data off the disk. 
#' 
#' A series of QC measures will be called when running this function. 
#' First, the rowData input must be either be a \code{GRanges} 
#' object or something than can be coerced into one (e.g. data frame or file).
#' The length of the rowData input after coercion must be equal to the maximum
#' index of the row in the sql data table. 
#' 
#' Next, the colData length must be equal to the maximum index of the column
#' in the sql data after coercion to a data.frame.
#' 
#' These maximums will be pulled on constuction, which ensures that the table
#' name, by default is 'data' exists in the sql file.
#' 
#' The sparse format suggests that there are three columns
#' in the sql table being referenced; two of which contain
#' 'row' and 'column' as labels that indicate the index of 
#' the row and column elements. 
#' 
#' @param rowData A \code{GRanges} object or something that can be coerced into one.
#' Type checking is performed in the call to \code{RangedSummarizedExperiment} 
#' @param colData A data frame of sample annotations
#' @param sqlfile String pointing to the sql file
#' @param lookupTableName = "data" Name of data table in sql file
#'
#' @return Returns a \code{sqlSparseHandle} object. 
#'
#' @examples
#' library(GenomicRanges)
#' nm <- c("chr", "start", "stop")
#' rt <- read.table(system.file("extdata", "dat1_row.bed", package = "yolo"))
#' rowData <- GRanges(setNames(rt, nm))
#' colData <- read.table(system.file("extdata", "dat1_col.txt", package = "yolo"))
#' sqlf <- system.file("extdata", "dat1.sqlite", package="yolo")
#' d1 <- sqlSparseHandleMake(rowData, colData, sqlfile=sqlf)
#' 
#' @import RSQLite
#' @import GenomicRanges
#' @import S4Vectors
#' @import SummarizedExperiment
#' @export
setGeneric(name = "sqlSparseHandleMake", def = function(rowData, colData, sqlfile, lookupTableName = "data")
    standardGeneric("sqlSparseHandleMake"))

#' @rdname sqlSparseHandleMake
setMethod("sqlSparseHandleMake", signature("ANY", "ANY", "character", "ANY"),
        definition = function(rowData, colData, sqlfile, lookupTableName = "data") {
        options(scipen=999)
            
        rowdim <- length(rowData)
        coldim <- dim(colData)[1]
            
        # Grab Max and Min; make sure that sqlfile exists
        con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=sqlfile)
        rowmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(row) from ', lookupTableName)))
        columnmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(column) from ', lookupTableName)))
        RSQLite::dbDisconnect(con)
            
        if(rowdim != rowmax) stop("Dimension of rows in sql file and rowData don't match!")
        if(coldim != columnmax) stop("Dimension of columns in sql file and colData don't match!")    
        
        # Marry file backend attributes
        colData$sqlFileHandle <- sqlfile
        colData$lookupTableName <- lookupTableName
        
        # Constructor
        rse <- SummarizedExperiment(rowRanges=rowData, colData=DataFrame(colData))
        rowMap <- seq_len(dim(rse)[1])
        names(rowMap) <- as.character(rowMap)
        colMap <- seq_len(dim(rse)[2])
        names(colMap) <- as.character(colMap)
        b <- new("sqlSparseHandle", rse, rowMap = rowMap, colMap = colMap)
        return(b)
})

