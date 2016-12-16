library(Matrix)

df1 <- read.table(system.file("extdata", "dat1.csv", package = "yolo"), sep = "," , header = TRUE)

a <- sparseMatrix(i = df1[,1], j = df1[,2], x = df1[,3])
sm <- summary(a)

# Compressed Column Format
colidx <- match(unique(sm$j), sm$j)
colmat <- sm[,c(1,3)]



# Compressed Row Format
sm <- sm[order(sm$i),]
rowidx <- match(unique(sm$i), sm$i)
rowmat <- sm[,c(2,3)]

rescue <- 1