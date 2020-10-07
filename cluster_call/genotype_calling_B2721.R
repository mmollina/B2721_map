#####
## Analytic procedures to construct the B2721 potato map - genotype calling
## -----------------------------------------------------------
## Author: Marcelo Mollinari
## Date: Sun Jul 19, 2020
## Bioinformatics Research Center
## North Carolina State University 
#####
## devtools::install_github("mmollina/ClusterCall")
## devtools::install_github(repo = "mmollina/mappoly", ref = "df2050e55787ff68c59ee900d6184bbbaf7d3996", dependencies = F)
require(ClusterCall)
require(mappoly)
dat.raw<-read.csv(file = "B2721_x_y.csv") ## Supplemenntary File S3
dat.raw[1:10,1:10]
dim(dat.raw)

# removing unused columns
datxy<-dat.raw[,-c(1,3,4)]

# spliting data in two matrices, one containing x information and other containing y
X<-grep(pattern = "X", colnames(datxy))
Y<-grep(pattern = "Y", colnames(datxy))
datx<-datxy[,X]
daty<-datxy[,Y]

# merging information from both parents
Px<-apply(datx[,1:2], 1, sum)
Py<-apply(daty[,1:2], 1, sum)
Qx<-apply(datx[,3:4], 1, sum)
Qy<-apply(daty[,3:4], 1, sum)

## Functions to transfrom the Cartesian coordinates to polar
raw2theta<-function(x,y)
  2*atan(y/x)/pi

raw2r<-function(x,y)
  sqrt(x^2 + y^2)

##Writting files containing `theta` and `r` information
## Theta
dat.theta<-matrix(NA, nrow(datx), ncol(datx)-2)
dat.theta[,1]<-raw2theta(Px, Py)
dat.theta[,2]<-raw2theta(Qx, Qy)
for(i in 5:ncol(datx))
{
  dat.theta[,i-2]<-raw2theta(datx[,i], daty[,i])
}
dimnames(dat.theta)<-list(as.character(datxy[,1]), c("Atlantic", "B1829.5", unlist(strsplit(colnames(datx), ".X"))[-c(1:4)]))
write.csv(x = dat.theta, file = "AtlanticxB1829_theta.csv")

## r
dat.r<-matrix(NA, nrow(datx), ncol(datx)-2)
dat.r[,1]<-raw2r(Px, Py)
dat.r[,2]<-raw2r(Qx, Qy)
for(i in 5:ncol(datx))
{
  dat.r[,i-2]<-raw2r(datx[,i], daty[,i])
}
dimnames(dat.r)<-list(as.character(datxy[,1]), c("Atlantic", "B1829.5", unlist(strsplit(colnames(datx), ".X"))[-c(1:4)]))
write.csv(x = dat.r, file = "AtlanticxB1829_r.csv")

## Running ClusterCall
ab <- read.pop(theta.file = "AtlanticxB1829_theta.csv", r.file = "AtlanticxB1829_r.csv")
AB <- CC.bipop(ab, parent1 = "Atlantic", parent2 = "B1829.5", n.core = 16)

#inspect.marker(AB, "solcap_snp_c1_10001")

## Genomic Info
## Updating genomic position
#### Get potato infinium 8303 sequences ####
solcap.snp.pos <- ape::read.gff(file = "potato_8303SNPs_potato_dm_v4.03.gff3")
head(solcap.snp.pos)
solcap.snp.pos$snp.names <- sapply(strsplit(solcap.snp.pos$attributes, split = ";|="), function(x) x[4])
solcap.snp.pos <- solcap.snp.pos[order(solcap.snp.pos$snp.names),]
solcap.snp.pos <- solcap.snp.pos[!duplicated(solcap.snp.pos$snp.names),]
solcap.snp.pos$ch <- sapply(strsplit(as.character(solcap.snp.pos$seqid), split = "chr"), function(x) as.numeric(x[2]))
solcap.snp.pos$ch[solcap.snp.pos$ch == 0] <- NA
solcap.snp.pos <- solcap.snp.pos[order(solcap.snp.pos[,"ch"],solcap.snp.pos[,"start"],decreasing=FALSE),]
solcap.snp.pos <- data.frame(chr = solcap.snp.pos$ch, pos = solcap.snp.pos$start, row.names = solcap.snp.pos$snp.names)
head(solcap.snp.pos)

## Filtering non conforming markers
pos<-ch<-geno<-nm<-nm.nc<-NULL
ct<-0
NC<-NULL
nc.limit <- 5
for(j in rownames(solcap.snp.pos))
{
  ct <- ct + 1
  cat(".")
  if(ct%%100==0) cat("\n")
  # if the markers has multiple entries in the ClusterCall result, ignore it 
  if(sum(AB@info$marker%in%j) != 1)
    next()
  new.name <- strsplit(j, "solcap_snp_")[[1]][2]
  parent.dose<-c(AB@geno[j,AB@parent1], AB@geno[j,AB@parent2])
  # if any of the parents have NA in their genotype for that marker, ignore it
  if(any(is.na(parent.dose)))
    next()
  # if the marker is monomorphyc, ignore it
  if(length(table(AB@geno[j,])) == 1)
    next()
  g<-which(segreg_poly(m = 4, AB@geno[j,AB@parent1], AB@geno[j,AB@parent2]) > 0) - 1
  nc<-!AB@geno[j, -c(1:2)]%in%g
  NC <- cbind(NC, nc)
  nm.nc <- c(nm.nc, j)
  ## if the number of non conforming individuals 
  ## is bigger than the threshold, discard the marker
  if(sum(nc) > ceiling(nc.limit * (ncol(ab@theta) - 2)/100))
    next()
  ## For-non conforming markers use NA
  if(sum(nc) != 0)
    AB@geno[j, c(FALSE, FALSE, nc)]<-NA
  nm<-c(nm, new.name)
  geno<-rbind(geno, AB@geno[j,])
  ch<-c(ch, solcap.snp.pos[j,"chr"])
  pos<-c(pos, solcap.snp.pos[j,"pos"])
}
NC<-t(NC)
dat<-data.frame(snp_name = nm,
           P1 = geno[,1], 
           P2 = geno[,2], 
           sequence = ch, 
           sequence_position = pos,
           geno[,-c(1:2)])
head(dat)
write.csv(x = dat, file = "B2721_CC.csv", row.names = FALSE)
save.image("all_cluster_call_data.rda")
sink("sessionInfo_cluster_call.txt")
sessionInfo()
sink()
