#####
## Analytical procedures to construct the B2721 potato map - meiotic configurations
## --------------------------------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
####
#### Loading required packages ####
require(mappoly)
require(tidyverse)
require(reshape2)
require(gtable)
require(grid)
require(igraph)
#### Loading functions and data ####
## setwd(dir = "meiosis/")
source("detect_meiotic_configuration_functions.R")
source("multiplot.R")
load(file = "../mapping/B2721_analysis.rda")
final.maps <- lapply(final.maps, filter_map_at_hmm_thres, thres.hmm = 0.00001)

#### Arguments ####
## Minimum homolog probability to consider the analyzed individual inherited it
hom.prob.thresh = 0.90
## Minimum segment length (in cM) to consider a crossing over
seg.length.thresh = 10
## Minimum percentage of information in a chromosome to consider it
perc.info = 60
## Minimum distance between two crossing overs to attribute exchange homolog pairs
dist.thresh = 2

#### Genotype probability computation: no step ####
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl,  "B2721")
genoprob.no.step <- parallel::parLapply(cl, final.maps, calc_genoprob_error, error = 0.05)
parallel::stopCluster(cl)

#### Recombination chains ####
## Gathering recombination points
RES <- vector("list", 12)
for(ch in 1:12){
  cat("   Chromosome",ch, "------------\n")
  gp <- genoprob.no.step[[ch]]
  all.ind <- dimnames(gp$probs)[[3]]
  c.p2<-c.p1<-res.p2<-g.p2<-o.p2<-res.p1<-g.p1<-o.p1<-vector("list", length(all.ind))
  names(c.p1)<-names(c.p2)<-names(res.p2)<-names(g.p2)<-names(o.p2)<-names(res.p1)<-names(g.p1)<-names(o.p1)<-all.ind
  for(individual in all.ind){
    cat("       ", individual, "\n")
    ## p1
    pr.hom <-  homo_prob(individual, gp, parent = 1)
    x1<-plot_recombination_points(pr.hom = pr.hom, 
                                  map = gp$map, 
                                  individual = individual, 
                                  parent = 1,
                                  hom.prob.thresh = hom.prob.thresh, 
                                  seg.length.thresh = seg.length.thresh, 
                                  perc.info = perc.info,
                                  map.mappoly = final.maps[[ch]], 
                                  dat.mappoly = B2721,
                                  dist.thresh = dist.thresh,
                                  title.plot =  individual)
    ## p1 Results
    res.p1[[individual]] <- x1$summary
    g.p1[[individual]] <- x1$plot
    o.p1[[individual]] <- x1$plot.orig
    c.p1[[individual]] <- x1$meiotic.configuration
    ##################
    ## p2
    pr.hom <-  homo_prob(individual, gp, parent = 2)
    x2<-plot_recombination_points(pr.hom = pr.hom, 
                                  map = gp$map, 
                                  individual = individual, 
                                  parent = 2,
                                  hom.prob.thresh = hom.prob.thresh, 
                                  seg.length.thresh = seg.length.thresh, 
                                  perc.info = perc.info,
                                  map.mappoly = final.maps[[ch]], 
                                  dat.mappoly = B2721,
                                  dist.thresh = dist.thresh,
                                  title.plot =  individual)
    ## p2 Results
    res.p2[[individual]] <- x2$summary
    g.p2[[individual]] <- x2$plot
    o.p2[[individual]] <- x2$plot.orig
    c.p2[[individual]] <- x2$meiotic.configuration
  }
  RES[[ch]] <- list(P1 = list(res.p1 = res.p1, g.p1 = g.p1, o.p1 = o.p1, c.p1 = c.p1), 
                    P2 = list(res.p2 = res.p2, g.p2 = g.p2, o.p2 = o.p2, c.p2 = c.p2))
}

## Getting recombination chain results
P1 <- sapply(RES, get_rec_chain, parent = 1)
round(100*sum(P1["Inconclusive",])/sum(P1),1)
## 7.3% were inconclusive
P2 <- sapply(RES, get_rec_chain, parent = 2)
round(100*sum(P2["Inconclusive",])/sum(P2),1)
## 9.8% were inconclusive

## Eliminating inconclusive cases. 
## Also summing 1 and 2 since in 
## both cases max there are two 
## homologs involved in the meiosis

## Parent 1
temp<-apply(P1[2:3,], 2, sum)
P1<-P1[-c(2,6), ]
P1[2,]<-temp
for(i in 1:12){
  P1[,i] <- round(100*P1[,i]/sum(P1[,i]),4)
}
## Parent 2
temp<-apply(P2[2:3,], 2, sum)
P2<-P2[-c(2,6),]
P2[2,]<-temp
for(i in 1:12){
  P2[,i] <- round(100*P2[,i]/sum(P2[,i]),4)
}
dimnames(P1) <- dimnames(P2) <- list(c("no c.o.","2","3","4"), paste0("ch", 1:12))
#dimnames(P1) <- dimnames(P2) <- list(c("no c.o.","1","2","3","4","Inc"), paste0("ch", 1:12))
df<-rbind(data.frame(parent = "Atlantic", reshape2::melt(P1)),
          data.frame(parent = "B1829", reshape2::melt(P2)))
df$valency <- ifelse(df$Var1 == "3" | df$Var1 == "4", "mul", "biv")  

png(filename = "~/repos/B2721_map/meiosis/meiotic_configurations.png", 
    width = 7.5, height = 5, units = "in", res = 200)
ggplot(data=df, aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") + facet_wrap( ~ parent, nrow = 2) +
  scale_fill_brewer(palette="Paired") +
  theme_minimal()
dev.off()

ggplot(data=df, aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") + facet_wrap( ~ parent, nrow = 2) +
  scale_fill_brewer(palette="Paired") +
  theme_minimal()
ggsave(filename = "meiotic_configurations.png", 
       width = 7.5, height = 5, units = "in", dpi = 200)
ggsave(filename = "meiotic_configurations.pdf", 
       width = 7.5, height = 5, units = "in")

#### Meiotic summary ####
round(P1,1)
round(P2,1)
P1.raw <- sapply(RES, get_rec_chain, parent = 1)
P2.raw <- sapply(RES, get_rec_chain, parent = 2)
#### Preferential Pairing ####
x1 <- calc_prefpair_profiles(genoprob)
png("pref_pairing.png", 
    width = 7.5, height = 7.5, units = "in", res = 200)
plot(x1, min.y.prof = 0.1, max.y.prof = .5, thresh = 0.01)
dev.off()

#### Population haplotype ####
h.prob.solcap<-calc_homoprob(genoprob)
print(h.prob.solcap)
plot(h.prob.solcap, ind = 1, lg = 1)
plot(h.prob.solcap, stack = TRUE, ind = 2, lg = 1:12)
head(h.prob.solcap$homoprob)
write.csv(h.prob.solcap$homoprob, file = "B2721_probabilistic_haplotypes.csv")
#### Saving image and session Info####
save(P1, P1.raw, P2, P2.raw, file = "meiosis_results.rda")
sink("sessionInfo_meiosis.txt")
sessionInfo()
sink()


