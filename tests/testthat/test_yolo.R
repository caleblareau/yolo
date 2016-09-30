
context(".rds of sqlSparseHandle are what we expect")
library(GenomicRanges)
nm <- c("chr", "start", "stop")

test_that("dat1.rds looks good", {
    rt <- read.table(system.file("extdata", "dat1_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat1_col.txt", package = "yolo"))
    d1 <- sqlSparseHandleMake(rowData, colData, sqlfile=system.file("extdata","dat1.sqlite",package="yolo"))
    #saveRDS(d1, file = "inst/rds/dat1_ssh.rds")
    drds <- readRDS(system.file("rds", "dat1_ssh.rds", package = "yolo"))
    expect_equal(dim(d1), dim(drds))
})

test_that("dat2.rds looks good", {
    rt <- read.table(system.file("extdata", "dat2_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat2_col.txt", package = "yolo"))
    d2 <- sqlSparseHandleMake(rowData, colData, sqlfile=system.file("extdata","dat2.sqlite",package="yolo"))
    #saveRDS(d2, file = "inst/rds/dat2_ssh.rds")
    drds <- readRDS(system.file("rds", "dat2_ssh.rds", package = "yolo"))
    expect_equal(dim(d2), dim(drds))
})

test_that("dat3.rds looks good", {
    rt <- read.table(system.file("extdata", "dat3_row.bed", package = "yolo"))
    rowData <- GRanges(setNames(rt, nm))
    colData <- read.table(system.file("extdata", "dat3_col.txt", package = "yolo"))
    d3 <- sqlSparseHandleMake(rowData, colData, sqlfile=system.file("extdata","dat3.sqlite",package="yolo"))
    #saveRDS(d3, file = "inst/rds/dat3_ssh.rds")
    drds <- readRDS(system.file("rds", "dat3_ssh.rds", package = "yolo"))
    expect_equal(dim(d3), dim(drds))
})