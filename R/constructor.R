#' @include classDefinitions.R
NULL

#' Create a \code{rseHandle} object
#'
#' \code{rseHandleMake} takes the a minimum of three parameters specified below
#' and creates a \code{rseHandle} object which can be subset and
#' manipulated without reading any data off the disk. At a minimum,
#' the rowData, colData, and backend file are needed. By default, we
#' assume that the file format will be an sqlite file and the format
#' will be sparse, but these can be differentially specified to an 
#' HDF5 file 
#' 
#' A series of QC measures will be called when running this function. 
#' First, the rowData input must be a \code{GRanges} object.
#' The length of the rowData input after coercion must be equal to the maximum
#' index of the row in the particular file/table. 
#' 
#' Next, the colData length must be equal to the maximum index of the column
#' in the linked file/table. 
#' 
#' These maximums will be pulled on constuction, which ensures that the table
#' name, by default is 'data' exists in the sql file.
#' 
#' The sparse format demands that there are three columns
#' in the file/table being referenced; two of which contain
#' 'row' and 'column' as labels that indicate the index of 
#' the row and column elements. 
#' 
#' The normal format demands that there are n columns in the
#' file/table being referenced where n is the numerb of samples
#' and the columns are named 
#' 
#' @param rowData A \code{GRanges} object or something that can be coerced into one.
#' Type checking is performed in the call to \code{RangedSummarizedExperiment} 
#' @param colData A data frame of per-sample annotations
#' @param lookupFileName String pointing to the file backend
#' @param lookupTableName = "data" Name of data table in file backend
#' @param lookupFileFormat = "sqlite" can also be "HDF5". The format of the
#' file providing the data backend of this constructed object. 
#' @param lookupFileType = "sparse" can also be "normal". The type of matrix
#' contained in the file/table. 
#'
#' @return Returns a \code{sqlSparseHandle} object. 
#'
#' @examples
#' library(GenomicRanges)
#' nm <- c("chr", "start", "stop")
#' rt <- read.table(system.file("extdata", "dat1_row.bed", package="yolo"))
#' rowData <- GRanges(setNames(rt, nm))
#' colData <- read.table(system.file("extdata", "dat1_col.txt", package="yolo"))
#' sqlf <- system.file("extdata", "dat1.sqlite", package="yolo")
#' d1 <- rseHandleMake(rowData, colData, lookupFileName=sqlf)
#' 
#' @import RSQLite
#' @import GenomicRanges
#' @import S4Vectors
#' @import SummarizedExperiment
#' @importFrom rhdf5 h5ls
#' @export
setGeneric(name = "rseHandleMake", def = function(rowData, colData, lookupFileName,
            lookupTableName = "data", lookupFileType = "sqlite", lookupFileFormat = "sparse")
    standardGeneric("rseHandleMake"))

#' @rdname rseHandleMake
setMethod("rseHandleMake", signature("ANY", "ANY", "character", "ANY"),
        definition = function(rowData, colData, lookupFileName, lookupTableName = "data",
                              lookupFileType = "sparse", lookupFileFormat = "sqlite") {
        options(scipen=999)
        
        stopifnot(lookupFileFormat %in% c("HDF5", "sqlite"))
        stopifnot(lookupFileType %in% c("normal", "sparse"))
        if(lookupFileType == "HDF5" & lookupFileFormat == "sparse") stop("HDF5/sparse combination not supported...")
            
        rowdim <- length(rowData)
        coldim <- dim(colData)[1]
            
        # Grab Max and Min ranges and ensure that file exists
        if(lookupFileFormat == "sqlite" & lookupFileType == "sparse"){
            con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=lookupFileName)
            rowmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(row) from ', lookupTableName)))
            columnmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(column) from ', lookupTableName)))
            dbDisconnect(con)
        } else if(lookupFileFormat == "HDF5" & lookupFileType == "normal"){
            rxc <- h5ls("dat.hdf5")[h5ls("dat.hdf5")$name == basename(lookupTableName), 5]
            if(length(rxc) != 1) stop("Couldn't discern a single table from specification; check HDF5 file")
            rc <- as.numeric(trimws(strsplit(rxc, "x")[[1]]))
            rowmax <- rc[1]
            columnmax <- rc[2]
        } else if(lookupFileFormat == "sqlite" & lookupFileType == "normal"){
            con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=lookupFileName)
            columnmax <- length(dbListFields(con, lookupTableName)) -1 
            sqlcmd <- paste0("SELECT Count(*) FROM ", lookupTableName)
            rowmax <- as.numeric(dbGetQuery(conn=con, statement=sqlcmd))
            dbDisconnect(con)
        } else {
            stop("Invalid specification of lookupFileFormat and lookupFileType")
        }
                
        if(rowdim != rowmax) stop("Dimension of rows in file/table and rowData don't match!")
        if(coldim != columnmax) stop("Dimension of columns in file/table and colData don't match!")    
        
        # Marry file backend attributes
        colData$lookupFileName   <- lookupFileName
        colData$lookupTableName  <- lookupTableName
        colData$lookupFileType   <- lookupFileType
        colData$lookupFileFormat <- lookupFileFormat
        
        # Constructor
        rse <- SummarizedExperiment(rowRanges=rowData, colData=DataFrame(colData))
        rowMap <- seq_len(dim(rse)[1])
        names(rowMap) <- as.character(rowMap)
        colMap <- seq_len(dim(rse)[2])
        names(colMap) <- as.character(colMap)
        b <- new("rseHandle", rse, rowMap = rowMap, colMap = colMap)
        return(b)
})

