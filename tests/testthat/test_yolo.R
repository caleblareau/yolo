# Sets of test for package

library(GenomicRanges)
library(rhdf5)

nm <- c("chr", "start", "stop")

context(".rds of Sparse SQlite Handles (ssh) are what we expect")

test_that("dat1.rds looks good", {
    rt <- read.table(system.file("extdata", "dat1_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat1_col.txt", package = "yolo"))
    d1 <- rseHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat1.sqlite",package="yolo"))
    #saveRDS(d1, file = "inst/rds/dat1_ssh.rds")
    dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    expect_equal(dim(d1), dim(dat1))
})

test_that("dat2.rds looks good", {
    rt <- read.table(system.file("extdata", "dat2_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat2_col.txt", package = "yolo"))
    d2 <- rseHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat2.sqlite",package="yolo"))
    #saveRDS(d2, file = "inst/rds/dat2_ssh.rds")
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    expect_equal(dim(d2), dim(dat2))
})

test_that("dat3.rds looks good", {
    rt <- read.table(system.file("extdata", "dat3_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat3_col.txt", package = "yolo"))
    d3 <- rseHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat3.sqlite",package="yolo"))
    #saveRDS(d3, file = "inst/rds/dat3_ssh.rds")
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    expect_equal(dim(d3), dim(dat3))
})

context("Subsetting works")

test_that("subsetting by range works in handles", {
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    reg <- GRanges(seqnames=c("chr1"),ranges=IRanges(start=c(3318000),end=c(3340000)))
    s1 <- subsetByOverlaps(as(dat3, "RangedSummarizedExperiment"), reg)
    s2 <- as(subsetByOverlaps(dat3, reg), "RangedSummarizedExperiment")
    expect_equal(s1, s2)
})

test_that("Square bracket operator behaves", {
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    d1 <- dat2[1:5, 3:4]
    d2 <- rseHandleSubset(dat2, 1:5, 3:4)
    expect_equal(d1, d2)
})

context("HDF5")
test_that("Adding two samples behaves", {
    dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    dat4 <- dat1 + dat3
    expect_equal(dim(dat4)[2], 35)
    expect_warning(expect_error(dat1 + dat2))
})


context("Addition")

test_that("Adding two samples behaves", {
    dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    dat4 <- dat1 + dat3
    expect_equal(dim(dat4)[2], 35)
    expect_warning(expect_error(dat1 + dat2))
})
