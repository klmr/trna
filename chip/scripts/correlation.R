source('../common/scripts/basic.R')
source('scripts/load-data.R')
local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/usage.R')
})

trnaGroupFamilyAndType <- function () {
    trnaByAcceptor <<- groupby(trnaNormDataCond, trnaAnnotation$Acceptor, sum)
    trnaByType <<- groupby(trnaNormDataCond, trnaAnnotation$Type, sum)
}
