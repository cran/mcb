library(mcb)

data(Diabetes) # load data
x <- Diabetes[,c('S1', 'S2', 'S3', 'S4', 'S5')]
y <- Diabetes[,c('Y')]
x <- data.matrix(x)
y <- data.matrix(y)

test_that("mcb works", {
  mcb(x=x, y=y)
})
