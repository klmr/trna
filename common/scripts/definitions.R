library(RColorBrewer)

source('../../../scripts/basic.R', chdir = TRUE)

tissues <- c(liver = 'liver', brain = 'brain')
stages <- c('e15.5', 'e18.5', 'P0.5', 'P4', 'P22', 'P29')
names(stages) <- stages

fullBrewer <- function (name) brewer.pal(brewer.pal.info[name, 'maxcolors'], name)

colors <- c("#597CCB", "#3D8E11", "#7C0D0C", "#C47E1F", "#603D71", "#CE4A92", "#A4B962")
colors <- c(colors, grey = '#4C4C4C')

contrastColors <- colorRampPalette(fullBrewer('PRGn'))(30)
# !!! The order of these colours is important, must reflect `tissues`.
tissueColor <- c(liver = colors[3], brain = colors[4])[tissues]

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
    sapply(str, function (str)
        paste(regswitch(unlist(strsplit(oldReadable(str), ' ')),
                              e15.5 = 'E15.5',
                              e18.5 = 'E18.5',
                              '[dD]o(\\d{4})' = 'DO\\1'),
              collapse = ' '))
