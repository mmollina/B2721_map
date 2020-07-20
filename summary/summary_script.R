#####
## Analitical procedures to construct the B2721 potato map - Map Summary
## -----------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
#####
require(mappoly)
load("~/repos/B2721_map/mapping/B2721_analysis.rda")
load("~/repos/B2721_map/meiosis/meiosis_results.rda")
setwd("summary/")

B2721
pval.bonf
100*length(mrks.chi.filt$exclude)/B2721$n.mrk
seq.unique

(S<-sum(c(seq.unique$seq.dose.p == 1 & seq.unique$seq.dose.q == 0,
          seq.unique$seq.dose.p == 1 & seq.unique$seq.dose.q == 4,
          seq.unique$seq.dose.p == 0 & seq.unique$seq.dose.q == 1,
          seq.unique$seq.dose.p == 4 & seq.unique$seq.dose.q == 1,
          seq.unique$seq.dose.p == 3 & seq.unique$seq.dose.q == 0,
          seq.unique$seq.dose.p == 3 & seq.unique$seq.dose.q == 4,
          seq.unique$seq.dose.p == 0 & seq.unique$seq.dose.q == 3,
          seq.unique$seq.dose.p == 4 & seq.unique$seq.dose.q == 3)))
(DS<-simplex<-sum(c(seq.unique$seq.dose.p == 1 & seq.unique$seq.dose.q == 1,
                    seq.unique$seq.dose.p == 1 & seq.unique$seq.dose.q == 3,
                    seq.unique$seq.dose.p == 3 & seq.unique$seq.dose.q == 3,
                    seq.unique$seq.dose.p == 3 & seq.unique$seq.dose.q == 1)))
(M<-simplex<-sum(c(seq.unique$seq.dose.p == 2 & seq.unique$seq.dose.q == 0,
                   seq.unique$seq.dose.p == 2 & seq.unique$seq.dose.q == 1,
                   seq.unique$seq.dose.p == 2 & seq.unique$seq.dose.q == 3,
                   seq.unique$seq.dose.p == 2 & seq.unique$seq.dose.q == 4,
                   seq.unique$seq.dose.p == 0 & seq.unique$seq.dose.q == 2,
                   seq.unique$seq.dose.p == 1 & seq.unique$seq.dose.q == 2,
                   seq.unique$seq.dose.p == 2 & seq.unique$seq.dose.q == 2,
                   seq.unique$seq.dose.p == 3 & seq.unique$seq.dose.q == 2,
                   seq.unique$seq.dose.p == 4 & seq.unique$seq.dose.q == 2)))
S;DS;M
round(100*S/(S + DS + M),1) ; round(100*DS/(S + DS + M),1) ; round(100*M/(S + DS + M),1)

png("filtered_data.png", 
    width = 7.5, height = 5, units = "in", res = 200)
plot(seq.unique)
dev.off()

all.rf.pairwise
round(100*sum(diag(grs$seq.vs.grouped.snp))/sum(grs$seq.vs.grouped.snp),1)

sink("map_summary.txt")
summary_maps(map.object = final.maps)
sink()

range(sapply(final.maps, function(x) round(sum(imf_h(x$maps[[1]]$seq.rf)),2)))
mean((sapply(final.maps, function(x) round(sum(imf_h(x$maps[[1]]$seq.rf)),2))))
png("final_maps.png", 
    width = 7.5, height = 7.5, units = "in", res = 200)
plot_map_list(final.maps, col = "ggstyle")
dev.off()

png("genome_vs_map.png", 
    width = 7.5, height = 7.5, units = "in", res = 200)
plot_genome_vs_map(final.maps, same.ch.lg = TRUE)
dev.off()
