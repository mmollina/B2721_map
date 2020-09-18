#####
## Analytic procedures to construct the B2721 potato map - Map Summary
## -----------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
#####
require(mappoly)
require(ggplot2)
load("../mapping/B2721_analysis.rda")
load("../meiosis/meiosis_results.rda")
source("../meiosis/detect_meiotic_configuration_functions.R")
setwd("summary/")
#### Initial Summary ####
## Data from ClusterCall
B2721
## Bonferroni's p.value
pval.bonf
## Percentage removed due to segregation distortion
100*length(mrks.chi.filt$exclude)/B2721$n.mrk

## Percentage removed due to segregetion distortion and redundancy  
round(100*(B2721$n.mrk-length(seq.unique$seq.num))/B2721$n.mrk,1)

## Number of markers per class
d<-table(seq.unique$seq.dose.p, seq.unique$seq.dose.q)
S <-  d["0","1"] + d["0","3"] + d["1","0"] + d["3","0"] +
      d["4","1"] + d["4","3"] + d["1","4"] + d["3","4"]
DS <- d["1","1"] + d["1","3"] + d["3","1"] + d["3","3"]
M <-  d["2","0"] + d["2","1"] + d["2","2"] + d["2","3"] + d["2","4"] +
      d["0","2"] + d["1","2"] + d["3","2"] + d["4","2"]
S;DS;M
round(100*S/(S + DS + M),2) ; round(100*DS/(S + DS + M),2) ; round(100*M/(S + DS + M),2)

#### SFigure-S2 ####
## Supplementary Figure S2 - Filtered data used to build the B2721 population map
png("SFigure-S2_filtered_data.png", 
    width = 7.5, height = 5, units = "in", res = 200)
plot(seq.unique)
dev.off()
## Pairwise recombination fraction
all.rf.pairwise
## Percentage of markers that coincide with respective reference chomosomes when groupping
## using the recombination fraction matrix
round(100*sum(diag(grs$seq.vs.grouped.snp))/sum(grs$seq.vs.grouped.snp),1)

#### STable-S1 ####
## Map Summary: Supplementary Table S1
sink("STable-S1_map_summary.txt")
summary_maps(map.object = final.maps)
sink()
## Range and mean
range(sapply(final.maps, function(x) round(sum(imf_h(x$maps[[1]]$seq.rf)),2)))
mean((sapply(final.maps, function(x) round(sum(imf_h(x$maps[[1]]$seq.rf)),2))))

#### SFigure-S3 ####
## Complete map figure: Supplementary Figure S3
png("SFigure-S3_complete_linkage_map.png", 
    width = 7.5, height = 7.5, units = "in", res = 200)
plot_map_list(final.maps, col = "ggstyle")
dev.off()
#### SFile-S5 ####
## Complete map file: Supplementary File S5
export_map_list(final.maps, file = "SFile-S5_complete_linkage_map.csv")
#### SFigure-S4 ####
## Genetic map versus reference genome": Supplementary Figure S4
png("SFigure-S4_map_vs_genome.png", 
    width = 7.5, height = 6.5, units = "in", res = 200)
plot_genome_vs_map(final.maps, same.ch.lg = TRUE)
dev.off()
#### SFigure-S5 ####
## Preferential pairing profiles: Supplementary Figure S5
x1 <- calc_prefpair_profiles(genoprob)
png("SFigure-S5_preferential_pairing.png", 
    width = 7.5, height = 5, units = "in", res = 200)
plot(x1, min.y.prof = 0.1, max.y.prof = .5, thresh = 0.01, P = "Atlantic", Q = "B1829-5")
dev.off()
#### Figure-2 ####
## Distribution of number of homologs involved in a recombination chain - Figure 2
df<-rbind(data.frame(parent = "Atlantic", reshape2::melt(P1)),
          data.frame(parent = "B1829-5", reshape2::melt(P2)))
df$valency <- ifelse(df$Var1 == "3" | df$Var1 == "4", "mul", "biv")  
ggplot(data=df, aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") + facet_wrap(parent ~ ., nrow = 2) +
  scale_fill_brewer(palette="Paired") +
  theme_minimal()
ggsave(filename = "meiotic_configurations.pdf", 
       width = 7.5, height = 5, units = "in")
## Edit names using Adobe Illustrator 

#### Meiosis Summary ####
round(100*sum(P1.raw["Inconclusive",])/sum(P1.raw),1)
## 7.3% were inconclusive
round(100*sum(P2.raw["Inconclusive",])/sum(P2.raw),1)
## 9.8% were inconclusive
##Remaining
P1[3, ] <- apply(P1[c(3,4),], 2, sum)
P1 <- P1[-4,]
rownames(P1)[3] <- "3-4" 
P2[3, ] <- apply(P2[c(3,4),], 2, sum)
P2 <- P2[-4,]
rownames(P2)[3] <- "3-4" 
cbind(P1, P2)
round(apply(cbind(P1, P2), 1, range),1)
round(apply(cbind(P1, P2), 1, mean),1)

