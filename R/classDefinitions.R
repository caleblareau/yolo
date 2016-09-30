#' @include yolo.R
NULL

#' sqlSparseHandle
#'
#' A class to represent sparse matrices stored in
#' and sql file where the column and row annotations
#' are wrapped in a \code{RangedSummarizedExperiment}
#' class. 
#'
#' @slot rowMap Named vector mapping indicies of
#' the current rows to indicies in the original file. 
#' @slot colMap Named vector mapping indicies of
#' the current columns to the indicies in the original file. 
#' 
#' @export
sqlSparseHandle <- setClass("sqlSparseHandle",
    contains="RangedSummarizedExperiment",
    representation = representation(
        rowMap="integer", 
        colMap="integer"
    )
)
