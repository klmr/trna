source('../common/scripts/basic.R')

# !!! IMPORTANT !!!
# The Ensembl Biomart archive server is extremely slow. Therefore this code is
# FOR EXPOSITION ONLY. Use the transcripts file shipped with the data for this
# code instead.
downloadTranscripts <- function (target) {
    require(biomaRt)
    ensMart <- useMart('ENSEMBL_MART_ENSEMBL',
                       'mmusculus_gene_ensembl',
                       'may2012.archive.ensembl.org')
    attributes <- c('ensembl_gene_id', 'ensembl_transcript_id', 'coding')
    transcripts <- getBM(attributes, mart = ensMart)
    fasta <- mapply(c,
                    do.call(sprintf, c('>%s|%s', as.list(transcripts[, 1 : 2]))),
                    splitLines(transcripts$coding, 60),
                    USE.NAMES = FALSE, SIMPLIFY = FALSE) %|% unlist
    writeLines(fasta, target)
}

splitLines <- function (str, lineLength, collapse = '\n')
    lapply(regmatches(str, gregexpr(sprintf('.{0,%d}', lineLength),
                                    str, perl = TRUE)),
           paste, collapse = '\n') %|% unlist

#downloadTranscripts('../common/data/Mus_musculus.NCBIM37.67.transcripts.fa')
