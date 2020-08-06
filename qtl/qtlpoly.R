#####
## Analytical procedures to map QTL using B2721 potato map
## -----------------------------------------------------------
## Author: Guilherme da Silva Pereira
## Date: Thu Aug 6, 2020
## Department of Genetics
## Luiz de Queiroz College of Agriculture
## University of São Paulo
#####

### Install and load packages
# install.packages("devtools")
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/varComp/varComp_0.2-0.tar.gz")
# devtools::install_version("sommer", version = "3.6", repos = "http://cran.us.r-project.org")
# devtools::install_github("guilherme-pereira/qtlpoly", upgrade = FALSE) 
library(qtlpoly)

### Prepare data
pheno <- read.table("pheno_asreml.txt", header = TRUE, row.names = 1); head(pheno)
# Coverting plant yield from pounds to kg
py <- which(substr(colnames(pheno), start = 1, stop = 2) == "PY")
pheno[,c(py)] <- pheno[,c(py)]/2.2046 
# Loading object with conditional probabilities
load("../mapping/B2721_analysis.rda") 
str(genoprob)
# Making sure genotype names match from pheno and geno objects
for(c in 1:length(genoprob)) dimnames(genoprob[[c]]$probs)[[3]] <- gsub(x = dimnames(genoprob[[c]]$probs)[[3]], pattern = ".", replacement = "-", fixed = TRUE)
# Reading genotypic and phenotypic data
data <- read_data(ploidy = 4, geno.prob = genoprob, pheno = pheno, step = 1)
print(data, detailed = FALSE)
save(data, file = "data.RData")

### Perform random-effect multiple interval mapping (REMIM)
## Score-based resampling method to assess genome-wide significance
data.sim <- simulate_qtl(data = data, mu = 0, h2.qtl = NULL, var.error = 1, n.sim = 1000, missing = TRUE, seed = 123) 
print(data.sim, detailed = TRUE)
save(data.sim, file="data_sim_null.RData")
score.null <- null_model(data = data.sim$results, n.clusters = 16, plot = NULL)
save(score.null, file="score_null.RData")

## Running REMIM analysis
remim.mod <- remim(data = data, w.size = 20, sig.fwd = 0.20, sig.bwd = 0.05, score.null = score.null, d.sint = 1.5, n.clusters = 16, plot = "remim_resampl")
print(remim.mod)
save("remim.mod", file="remim_mod_resampl.RData")

## Fitting final models with detected QTL
fitted.mod <- fit_model(data = data, model = remim.mod, probs = "joint", polygenes = "none")
summary(fitted.mod)
save("fitted.mod", file="fitted_mod_resampl.RData")

## Estimating allele effects
est.effects <- qtl_effects(ploidy = 4, fitted = fitted.mod)
save("est.effects", file="est_effects_resampl.RData")

### Perform fixed-effect interval mapping (FEIM)
## Permutation method to assess genome-wide significance
perm <- permutations(data = data, n.sim = 1000, n.clusters = 16)
sig.lod <- numeric(length(perm$results)); for(p in 1:length(perm$results)) sig.lod[p] <- quantile(sort(perm$results[[p]]), 0.95)

## Running FEIM analysis
feim.mod <- feim(data = data, w.size = 20, sig.lod = sig.lod)
print(feim.mod)
save("perm","feim.mod", file="feim.RData")

sink("sessionInfo_mapping.txt")
sessionInfo()
sink()


