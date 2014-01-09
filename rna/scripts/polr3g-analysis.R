source('../common/scripts/basic.R')
source('scripts/load-data.R')

mrnaLoadData()
mrnaNormalizeData()

idx <- grep('Polr3g', mrnaAnnotation$Name)
pol3g <- cbind(data.frame(Gene = mrnaAnnotation[idx, 'Name']),
               mrnaNormDataCond[idx, ])

pol3g <- cbind(Gene = c(as.character(pol3g$Gene), 'Sum'),
               rbind(pol3g[-1], colSums(pol3g[-1])))

ymax <- round(max(pol3g[-1]), -2)
colors <- rainbow(3)
legendpos <- c(liver = 'bottomleft', brain = 'topleft')

for (tissue in tissues) {
    plot(indices(stages), type = 'n', ylim = c(0, ymax),
         xlab = 'Stages', ylab = 'Expression', xaxt = 'n',
         main = paste('Pol3g expression over time in', tissue))
    axis(1, at = indices(stages), labels = stages)
    for (i in 1 : 3)
        lines(indices(stages), pol3g[i, grep(tissue, colnames(pol3g))],
              col = colors[i])
    legend(legendpos[tissue], legend = pol3g$Gene, fill = colors, border = FALSE, bty = 'n')
}
