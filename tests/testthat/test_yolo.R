# Sets of test for package

library(GenomicRanges)
library(rhdf5)

nm <- c("chr", "start", "stop")

context(".rds of Sparse SQlite rseHandles (ssr) are what we expect")

test_that("dat1.rds looks good", {
    rt <- read.table(system.file("extdata", "dat1_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat1_col.txt", package = "yolo"))
    d1 <- yoloHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat1.sqlite",package="yolo"))
    se1 <-  yoloHandleMake(rt, colData, lookupFileName=system.file("extdata","dat1.sqlite",package="yolo"))
    #saveRDS(d1, file = "inst/rds/dat1_ssh.rds")
    dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    expect_equal(dim(d1), dim(dat1))
})

test_that("dat2.rds looks good", {
    rt <- read.table(system.file("extdata", "dat2_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat2_col.txt", package = "yolo"))
    d2 <- yoloHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat2.sqlite",package="yolo"))
    #saveRDS(d2, file = "inst/rds/dat2_ssh.rds")
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    expect_equal(dim(d2), dim(dat2))
})

test_that("dat3.rds looks good", {
    rt <- read.table(system.file("extdata", "dat3_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat3_col.txt", package = "yolo"))
    d3 <- yoloHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat3.sqlite",package="yolo"))
    #saveRDS(d3, file = "inst/rds/dat3_ssh.rds")
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    expect_equal(dim(d3), dim(dat3))
})


context(".rds of Sparse SQlite seHandles are what we expect")
test_that("se1.rds looks good", {
    rt <- read.table(system.file("extdata", "dat1_row.bed", package = "yolo"))
    rowData <- setNames(rt, nm)
    colData <- read.table(system.file("extdata", "dat1_col.txt", package = "yolo"))
    se1 <-  yoloHandleMake(rt, colData, lookupFileName=system.file("extdata","dat1.sqlite",package="yolo"))
    #saveRDS(se1, file = "inst/rds/se1.rds")
    s1 <- readRDS(system.file("rds", "se1.rds", package = "yolo"))
    expect_equal(dim(s1), dim(se1))
    expect_equal(as.character(class(s1)), "seHandle")
})

test_that("se2.rds looks good", {
    rt <- read.table(system.file("extdata", "dat2_row.bed", package = "yolo"))
    rowData <- setNames(rt, nm)
    colData <- read.table(system.file("extdata", "dat2_col.txt", package = "yolo"))
    se2 <- yoloHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat2.sqlite",package="yolo"))
    #saveRDS(se2, file = "inst/rds/se2.rds")
    s2 <- readRDS(system.file("rds", "se2.rds", package = "yolo"))
    expect_equal(dim(s2), dim(se2))
    expect_equal(as.character(class(s2)), "seHandle")
})

test_that("se3.rds looks good", {
    rt <- read.table(system.file("extdata", "dat3_row.bed", package = "yolo"))
    rowData <- setNames(rt, nm)
    colData <- read.table(system.file("extdata", "dat3_col.txt", package = "yolo"))
    se3 <- yoloHandleMake(rowData, colData, lookupFileName=system.file("extdata","dat3.sqlite",package="yolo"))
    #saveRDS(se3, file = "inst/rds/se3.rds")
    s3 <- readRDS(system.file("rds", "se3.rds", package = "yolo"))
    expect_equal(dim(s3), dim(se3))
    expect_equal(as.character(class(s3)), "seHandle")
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
    d2 <- yoloHandleSubset(dat2, 1:5, 3:4)
    expect_equal(d1, d2)
})


context("Addition")

test_that("Adding two rse samples behaves", {
    dat1 <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    dat2 <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    dat3 <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    dat4 <- dat1 + dat3
    expect_equal(dim(dat4)[2], 35)
    expect_warning(expect_error(dat1 + dat2))
})

test_that("Adding two se samples behaves", {
    se1 <- readRDS(system.file("rds", "se1.rds", package = "yolo"))
    se2 <- readRDS(system.file("rds", "se2.rds", package = "yolo"))
    se3 <- readRDS(system.file("rds", "se3.rds", package = "yolo"))
    se4 <- se1 + se3
    expect_equal(dim(se4)[2], 35)
    expect_error(se1 + se2)
})