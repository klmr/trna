library(RColorBrewer)

source('../../../scripts/basic.R', chdir = TRUE)

tissues <- c(liver = 'liver', brain = 'brain')
stages <- c('e15.5', 'e18.5', 'P0.5', 'P4', 'P22', 'P29')
names(stages) <- stages

colors <- brewer.pal(8, 'Dark2')
# !!! The order of these colours is important, must reflect `tissues`.
tissueColor <- c(liver = '#D01B24', brain = '#E6AB02')

plotFamily <- 'Helvetica'

source('plot-matrix.R')
source('plot-pca.R')
source('plot-pairwise.R')
source('plot-correlations.R')

loadAminoAcids <- function () {
    aminoAcidPath <- '../common/data/amino_acids.tsv'
    aminoAcids <<- read.table(aminoAcidPath,
                              col.names = c('Long', 'Short'),
                              stringsAsFactors = FALSE)
}

loadGeneticCode <- function () {
    geneticCodeFile <- '../common/data/genetic_code.tsv'
    geneticCode <<- read.table(geneticCodeFile, row.names = 1,
                               col.names = c('', 'AA'),
                               stringsAsFactors = FALSE)
}

# FIXME Remove the TeXy hack once we switch to modules.
oldReadable <- readable
readable <- function (str)
    paste(regswitch(unlist(strsplit(oldReadable(str), ' ')),
                          e15.5 = 'E15.5',
                          e18.5 = 'E18.5',
                          '[dD]o(\\d{4})' = 'DO\\1'),
          collapse = ' ')
