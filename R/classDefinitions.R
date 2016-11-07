#' @include yolo.R
NULL

#' seHandle
#'
#' A class to represent  matrices stored in HDF5
#' and sqlite files where the column and row annotations
#' are wrapped in a \code{SummarizedExperiment}
#' class. 
#'
#' @slot rowMap Named vector mapping indicies of
#' the current rows to indicies in the original file. 
#' @slot colMap Named vector mapping indicies of
#' the current columns to the indicies in the original file. 
#' 
#' @export
seHandle <- setClass("seHandle",
    contains="SummarizedExperiment",
    representation = representation(
        rowMap="integer", 
        colMap="integer"
    )
)

#' rseHandle
#'
#' A class to represent  matrices stored in HDF5
#' and sqlite files where the column and row annotations
#' are wrapped in a \code{RangedSummarizedExperiment}
#' class. 
#'
#' @slot rowMap Named vector mapping indicies of
#' the current rows to indicies in the original file. 
#' @slot colMap Named vector mapping indicies of
#' the current columns to the indicies in the original file. 
#' 
#' @export
rseHandle <- setClass("rseHandle",
    contains="RangedSummarizedExperiment",
    representation = representation(
        rowMap="integer", 
        colMap="integer"
    )
)

setClassUnion("yoloHandle", c("seHandle", "rseHandle"))
