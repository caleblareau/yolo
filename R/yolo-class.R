#' @include sparseHiC.R
NULL

check_hic <- function(object) {
    errors <- character()
    names <- names(object@resolutionNamedList)
    if (!is.character(names)) {
        msg <- "Names of resolutions should be a character vector"
        errors <- c(errors, msg)
    }
    if (!is.character(object@sampleName)) {
        msg <- "Sample name should be a character"
        errors <- c(errors, msg)
    }
    if(any(is.na(suppressWarnings(as.numeric(names))))){
        msg <-"At least one resolution couldn't be coerced into a number"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}


#' A class to Hi-C data from a single sample but multiple resolutions
#'
#' @slot sampleName The name of the single sample that is the feature
#' of the sparseHiCdatum.
#' @slot resolutionNamedList A list of lists of Hi-C matrices named using the
#' string of the resolution.
#' @slot metaData A slot containing data corresponding to annotation data for
#' the particular sample.
#' @export
sparseHiCdatum <- setClass("sparseHiCdatum",
                           slots = c(sampleName = "character",
                                     resolutionNamedList = "list",
                                     metaData = "data.frame"),
                           validity = check_hic)


checkHiCdtm <- function(m) class(m)[1] == "sparseHiCdatum"

check_hic_samples <- function(object) {
    errors <- character()
    names <- names(object@HiCSamplesList)
    if (!is.character(names)) {
        msg <- "Names of samples should be a character vector"
        errors <- c(errors, msg)
    }
    #if(!all(sapply(object@HiCSamplesList, checkHiCdtm))){
    #    msg <- "Attempting to assemble one or more non-sparseHiCdatum objects"
    #    errors <- c(errors, msg)
    #}
    if (length(errors) == 0) TRUE else errors
}


#' A class to Hi-C data from multiple samples; collection of the sparseHiCdatum
#' objects
#'
#' @slot HiCSamplesList A list (named by samples) of lists (per sample resolutions)
#' of lists (Hi-C matrices)
#' @slot metaData A data.frame of data that may be useful for 
#' 
#' @export
sparseHiCdata <- setClass("sparseHiCdata",
                           slots = c(HiCSamplesList = "list", 
                                     metaData = "data.frame"),
                           validity = check_hic_samples)


