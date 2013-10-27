source('../common/scripts/basic.R')
source('scripts/load-data.R')
local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/correlation.R')
})

trnaGroupFamilyAndType <- function () {
    trnaByAcceptor <<- groupby(trnaNormDataCond, trnaAnnotation$Acceptor, mean)
    trnaByType <<- groupby(trnaNormDataCond, trnaAnnotation$Type, mean)
}

relativeData <- function (data) {
    sums <- apply(data, COLS, sum)
    t(apply(data, ROWS, function (row) row / sums))
}
