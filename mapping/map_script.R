#####
## Analitical procedures to construct the B2721 potato map
## -----------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
#####
## NOTE:to reproduce exactly the results presented in Pereira et al. (2020), install MAPpoly from the following Git Commit
## devtools::install_github(repo = "mmollina/mappoly", ref = "affa31d531e5c460386277934634e9e3c26b4e90", dependencies = F)
#### Functions ####
phasing_and_hmm_rf<-function(X){
  fl<-paste0("output_map_ch_", X$seq$sequence[1], ".txt")
  sink(fl)
  map<-est_rf_hmm_sequential(input.seq = X$seq,
                             start.set = 5,
                             thres.twopt = 10,
                             thres.hmm = 10,
                             extend.tail = 200,
                             twopt = X$tpt,
                             verbose = TRUE,
                             tol = 10e-3,
                             tol.final = 10e-4,
                             phase.number.limit = 60,
                             sub.map.size.diff.limit = 8,
                             info.tail = TRUE,
                             reestimate.single.ph.configuration = TRUE) 
  sink()
  return(map)
}
error_model<-function(X, tol = 10e-4){
  X$maps[[1]]$seq.rf <- rep(0.01, length(X$maps[[1]]$seq.rf))
  x<-est_full_hmm_with_global_error(input.map = X, 
                                    error = 0.05, 
                                    tol = tol, 
                                    verbose = FALSE)
  return(x)
}
require(mappoly)
require(ggplot2)
#### Initial analysis ####
## Loading dataset from ClusterCall
B2721 <- read_geno_csv(file.in  = "../cluster_call/B2721_CC.csv", ploidy = 4, elim.redundant = FALSE)
B2721
## Segregation test
pval.bonf <- 0.05/B2721$n.mrk
mrks.chi.filt <- filter_segregation(B2721, 
                                    chisq.pval.thres =  pval.bonf, 
                                    inter = FALSE)
seq.init<-make_seq_mappoly(mrks.chi.filt)

## Elinating redundant markers
seq.redundant <- elim_redundant(input.seq = seq.init)
print(seq.redundant)
plot(seq.redundant)
seq.unique <- make_seq_mappoly(seq.redundant)

png("figures_during_analysis/filtered_data.png", 
    width = 7.5, height = 5, units = "in", res = 200)
plot(seq.unique)
dev.off()

#### Two-point analysis ####
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.unique, 
                                   n.clusters = 24)
save(all.rf.pairwise, file = "twopt.RData")
#load("twopt.RData")

## Recombinarion fraction matrix
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
png("figures_during_analysis/full_rf_mat.png", 
    width = 7.5, height = 7.5, units = "in", res = 200)
plot(mat, ord = rownames(get_genomic_order(seq.unique)), index = FALSE)
dev.off()

#### Groupping ####
grs <- group_mappoly(input.mat = mat,
                     expected.groups = 12,
                     comp.mat = TRUE, 
                     inter = TRUE)
grs

## Correspondence with genome
z<-as.numeric(colnames(grs$seq.vs.grouped.snp)[1:12])

#### Assebling linkage groups (order based on genome) ####
LGS<-vector("list", 12)
for(ch in 1:12){
  cat("\n ~~~~~~ ch:", ch, "...\n")
  lg <- which(z==ch)
  s<-make_seq_mappoly(grs, lg)
  s<-make_seq_mappoly(B2721, names(which(s$sequence == ch)))
  tpt <- make_pairs_mappoly(all.rf.pairwise, input.seq = s)
  LGS[[ch]] <- list(seq = s,
                    tpt = tpt)
}
#### Parallel map construction and genotype probability computation ####
cl <- parallel::makeCluster(12)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterExport(cl,  "B2721")
MAPs.geno <- parallel::parLapply(cl, LGS, phasing_and_hmm_rf)
final.maps <- parallel::parLapply(cl, MAPs.geno, error_model)
genoprob <- parallel::parLapply(cl, final.maps, calc_genoprob_error, step = 1, error = 0.05)
parallel::stopCluster(cl)
#### Saving image and session Info####
save.image("B2721_analysis.rda")
sink("sessionInfo_mapping.txt")
sessionInfo()
sink()

