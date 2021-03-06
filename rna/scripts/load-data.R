mrnaDataFile <- '../common/data/genes.raw.htseq2.tsv'
mrnaMappingFile <- '../common/data/rna-seq-summary.tsv'
mrnaAnnotationFile <- '../common/data/gene.annot.tsv'

mrnaLoadData <- function () {
    if (exists('mrnaRawCounts'))
        return()

    mrnaRawCounts <- read.table(mrnaDataFile, header = TRUE, row.names = 1)
    mrnaMapping <- read.table(mrnaMappingFile, header = FALSE, row.names = 3,
                             col.names = c('Tissue', 'Stage', ''))
    mrnaMapping$Condition <- paste(mrnaMapping$Tissue, mrnaMapping$Stage, sep = '-')
    # Ensure corresponding order.
    mrnaMapping <- mrnaMapping[colnames(mrnaRawCounts), ]

    mrnaAnnotation <- read.table(mrnaAnnotationFile, header = TRUE)
    rownames(mrnaAnnotation) <- mrnaAnnotation$ID
    mrnaAnnotation$ID <- NULL
    mrnaAnnotation$GO <- NULL
    mrnaAnnotation$GOterm <- NULL
    mrnaAnnotation$KEGG <- NULL
    nr <- nrow(mrnaAnnotation)
    nonMitochondrial <- ! boolmask(grep('^MT:\\d+..\\d+$', mrnaAnnotation$locus), nr)
    nonFragment <- ! boolmask(grep('^NT_\\d+:\\d+..\\d+$', mrnaAnnotation$locus), nr)
    nonSex <- ! boolmask(grep('^[XY]:', mrnaAnnotation$locus), nr)
    proteinCoding <- grepl('protein_coding', mrnaAnnotation$source, fixed = TRUE)
    mrnaAnnotation <- mrnaAnnotation[nonMitochondrial & nonFragment & nonSex & proteinCoding, ]
    mrnaAnnotation$source <- NULL
    mrnaAnnotation$Chr <- sub(':.*', '', mrnaAnnotation$locus)
    mrnaAnnotation$Start <- as.numeric(sub('.*:(\\d+)\\..*', '\\1', mrnaAnnotation$locus))
    mrnaAnnotation$End <- as.numeric(sub('.*\\.(\\d+)', '\\1', mrnaAnnotation$locus))
    mrnaAnnotation$locus <- NULL

    mrnaRawCounts <- mrnaRawCounts[rownames(mrnaAnnotation), ]

    mrnaRawCounts <<- mrnaRawCounts
    mrnaMapping <<- mrnaMapping
    mrnaAnnotation <<- mrnaAnnotation
}

mrnaNormalizeData <- function () {
    if (exists('mrnaNormData'))
        return()

    require(DESeq)
    mrnaCds <- newCountDataSet(mrnaRawCounts, mrnaMapping$Condition)
    mrnaCds <- estimateSizeFactors(mrnaCds)
    mrnaNormData <- as.data.frame(counts(mrnaCds, normalized = TRUE))
    conditions <- unique(mrnaMapping$Condition)

    mrnaNormDataCond <- matrix(nrow = nrow(mrnaNormData), ncol = length(conditions))
    colnames(mrnaNormDataCond) <- conditions
    rownames(mrnaNormDataCond) <- rownames(mrnaNormData)

    # Same as `t(groupby(t(mrnaNormData), mrnaMapping$Condition, mean))` but
    # much faster. I need a transposed version of that function.
    for (cond in conditions) {
        replicates <- mrnaNormData[, rownames(mrnaMapping[mrnaMapping$Condition == cond, ])]
        if (! is.data.frame(replicates))
            mrnaNormDataCond[, cond] <- replicates
        else
            mrnaNormDataCond[, cond] <- rowMeans(replicates)
    }

    mrnaNormData <<- mrnaNormData
    mrnaNormDataCond <<- as.data.frame(mrnaNormDataCond)
}

mkdir('results')
mkdir('plots')
