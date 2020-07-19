## Declaring functions
#rm(list = ls())
#load("~/repos/Potato_Atlantic_B1829/mapping/genoprob_data_and_map_info.RData")
#genoprob <- calc_genoprob_error(input.map = final.maps[[1]], error = 0.05, verbose = TRUE)

## ind <- "B2721.026"
## parent <- 2
## map <- genoprob$map
## pr.hom <-  homo_prob(ind, genoprob, parent)
## image(pr.hom)
homo_prob<-function(ind, genoprob, parent){
  v <- c(4,36,400,4900)
  names(v) <- c(2,4,6,8)
  m <- as.numeric(names(which(dim(genoprob$probs)[1] == v)))
  ## genotype names
  geno.names<-dimnames(genoprob$probs)[[1]]
  ## conditional probability matrix
  pr.mat<-genoprob$probs[,,ind]
  ## choosing parent
  if(parent == 1 ){
    lt<-letters[1:m]
  } else if(parent == 2){
    lt<-letters[(m+1):(2*m)]
  } else stop("parents should be '1' or '2'.")
  #####
  ## probability matrix for the Homologs chromosomes
  pr.hom<-NULL
  for(i in 1:m)
  {
    y<-apply(pr.mat[grep(pattern = lt[i], geno.names),], 2, sum)
    pr.hom<-rbind(pr.hom, y)
  }
  rownames(pr.hom)<-lt
  return(pr.hom)  
}

## plot_homo_prob(pr.hom, map, ind)
plot_homo_prob<-function(pr.hom, map, ind, title.plot = NULL){
  if("a"%in%dimnames(pr.hom)[[1]]){
    parent <- 1    
  } else parent <- 2
  m <- dim(pr.hom)[1]
  if(is.null(title.plot))
    title.plot <- paste("Parent", parent)
  ## Transforming probability matrix for the Homologs chromosomes into a dataframe
  lt<-rownames(pr.hom)
  if(parent == 1){
    pal<-c("#800000", "#9A6324", "#808000", "#e6194B", "#f58231", "#ffe119", '#bcf60c', "#3cb44b")
    pal<-pal[ceiling(seq(1, 8, length.out = m))]
  } else {
    pal<-c('#aaffc3',"#469990", "#42d4f4", "#4363d8", "#000075", "#911eb4", '#e6beff', '#f032e6')
    pal<-pal[ceiling(seq(1, 8, length.out = m))]
  }
  df.pr<-NULL
  for(i in 1:m)
    df.pr<-rbind(df.pr, data.frame(pr = pr.hom[i,], 
                                   hom = lt[i], 
                                   position = map[names(pr.hom[i,])]))
  ## ggplot graphic
  p <- ggplot(df.pr, aes(x = position, y = pr, fill = hom, color  = hom)) +
    geom_density(stat = "identity", alpha = 0.7) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
    ggtitle(title.plot) + facet_grid(rows = vars(hom)) + theme_minimal() + ylab(label = "Homologs probabilty")
  return(p)
}

## hom.prob.thresh <- 0.8
## bool.pr.hom <- elim_low_prob_homo_segments(pr.hom, hom.prob.thresh)
## image(bool.pr.hom)
elim_low_prob_homo_segments<-function(pr.hom, hom.prob.thresh){
  m <- dim(pr.hom)[1]
  out<-apply(pr.hom, 2, function(x) {
    y<-rep(0,m)
    y[x > hom.prob.thresh]<-1
    y})
  rownames(out)<-rownames(pr.hom)
  # At least one Homologs has pr > hom.prob.thresh
  out <- out[,(apply(out, 2, sum) >= 1)]
  return(out)
}

## segment.lengths <- get_segment_length(bool.pr.hom, map, is.gap = FALSE)
## segment.lengths
get_segment_length<-function(bool.pr.hom, map, is.gap = FALSE){
  lt<-rownames(bool.pr.hom)
  A <- bool.pr.hom
  if(!is.gap){
    A[A == 1] <- 2
    A[A == 0] <- 1
    A[A == 2] <- 0
  }
  dimnames(A)[[1]]<-lt
  B0<-apply(A, 1, function(x) abs(diff(x)))
  D0<-apply(B0, 1, function(x) which(x==1))
  D0 <- B0[which(sapply(D0, length) > 0), ,drop=FALSE]
  B1<-vector("list", length(lt))
  names(B1)<-lt
  for(j in 1:length(lt)){
    x <- c(1, which(diff(A[j,])!=0), ncol(A))
    id<-A[j,x[-length(x)]+1]==0
    x <- cbind(x[-length(x)], x[-1])
    xlen <- NULL
    for(i in 1:nrow(x))
      xlen <- c(xlen, map[x[i,2]] - map[x[i,1]])
    B1[[j]] <- cbind(x, xlen)[id,]
    B1[[j]]<-matrix(B1[[j]], ncol = 3, byrow = FALSE)
    colnames(B1[[j]])<-c("start", "end", "length")
  }
  n.breaks<-matrix(apply(B0, 2, sum), nrow = 1, dimnames = list("ind", lt))
  list(segment.lengths = B1, n.breaks = n.breaks)
}

## seg.length.thresh <- 10
## bool.pr.hom.filt <- elim_small_and_incomplete_segments(bool.pr.hom, segment.lengths, seg.length.thresh)
## image(bool.pr.hom.filt)
elim_small_and_incomplete_segments<-function(bool.pr.hom, segment.lengths, seg.length.thresh){
  lt<-rownames(bool.pr.hom)
  for(j in 1:length(lt)){
    for(i in 1:nrow(segment.lengths$segment.lengths[[j]]))
      if(nrow(segment.lengths$segment.lengths[[j]]) > 0)
        if(segment.lengths$segment.lengths[[j]][i,3] < seg.length.thresh)
          bool.pr.hom[,segment.lengths$segment.lengths[[j]][i,1]:segment.lengths$segment.lengths[[j]][i,2]]<-NA
  }
  idrem<-which(apply(bool.pr.hom, 2, function(x) all(is.na(x))))
  if(length(idrem) != 0){
    bool.pr.hom.NA <- bool.pr.hom
    bool.pr.hom<- bool.pr.hom[ ,-idrem]
  } 
  mrk.names<-bool.pr.hom.res<-NULL
  for(i in 1:ncol(bool.pr.hom)){
    if(!all(is.na(bool.pr.hom[,i])))
      # All three Homologs have to have high probabilities, otherwise, we estimate using HMM + high quality markers
      if(sum(bool.pr.hom[,i], na.rm = TRUE) == length(lt)/2){
        bool.pr.hom.res <- cbind(bool.pr.hom.res, bool.pr.hom[,i])
        mrk.names<-c(mrk.names, colnames(bool.pr.hom)[i])
      }
  }
  colnames(bool.pr.hom.res)<-mrk.names
  return(bool.pr.hom.res)
}

## image(pr.hom)
## image(bool.pr.hom)
## image(bool.pr.hom.filt)
## co.pos1 <- detect_co(pr.hom, map)
## co.pos2 <- detect_co(bool.pr.hom, map)
## co.pos3 <- detect_co(bool.pr.hom.filt, map)
detect_co<-function(pr.hom, map, dist.thresh = 5){
  if(is.null(pr.hom))
    return(list(result = "Inconclusive", result.par = "Inconclusive"))
  lt<-rownames(pr.hom)
  ## Detecting crossing overs
  A<-apply(pr.hom, 2, function(x) {
    y<-rep(0,length(lt))
    y[order(x)[1:(length(lt)/2)]]<-1
    y})
  dimnames(A)[[1]]<-lt
  B<-apply(A, 1, function(x) abs(diff(x)))
  D<-apply(B, 1, function(x) which(x==1))
  D <- B[which(sapply(D, length) > 0), ,drop=FALSE]
  result<-"Zero"
  if(nrow(D)!=0){
    hom.pair <- apply(D, 1, function(x) paste0(lt[which(x==1)], collapse = ""))
    position <- round(map[names(hom.pair)],1) 
    z1 <- apply(D, 1, function(x) {y <- which(x==1); if(length(y)==2) return(y) else return(c(NA, NA))})
    from <- z1[1,]
    to <- z1[2,]
    point.id<-match(names(hom.pair), names(map))
    result<-data.frame(hom.pair, position, from, to, point.id)
  }
  if(as.character(result)[1]=="Zero")
    return(list(result = "Zero", result.par = "Zero"))
  id<-is.na(result$from)
  rem = NULL
  if(any(id)){
    rem<-result[which(id),]
    result<-result[-which(id),]
  }
  if(nrow(result)==0)
    return(list(result = "Inconclusive", result.par = "Inconclusive"))
  rem <- which(diff(result$position) < dist.thresh)
  if(length(rem) > 0){
    for(i in 1:length(rem))
      rem <- c(rem, rem[i] + 1)
    result <- result[-rem,,drop = FALSE]
    if(nrow(result)==0)
      return(list(result = "Inconclusive", result.par = "Inconclusive"))
  }
  return(list(result = result, result.par = result[,1:2]))
}


## gen.new <- calc_genoprob_err_one_ind(input.map  = map.mappoly, dat.err = dat.mappoly,  individual = individual) 
calc_genoprob_err_one_ind<-function(input.map, dat.err, individual, step = 0, phase.config = "best", error = 0.01, 
                                    th.prob = 0.95, restricted = TRUE, verbose = TRUE) 
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  m<-input.map$info$m
  dat.err$n.ind <- 1
  dat.err$ind.names <- individual
  dat.err$geno.dose <- dat.err$geno.dose[, individual, drop = FALSE]
  n.ind<-dat.err$n.ind
  Dtemp<-dat.err$geno.dose[input.map$info$mrk.names,,drop = FALSE]
  map.pseudo <- create_map(input.map)
  mrknames<-names(map.pseudo)
  n.mrk<-length(map.pseudo)
  indnames<-colnames(Dtemp)
  D<-matrix(m+1, nrow = length(map.pseudo), ncol = ncol(Dtemp), 
            dimnames = list(mrknames, indnames))
  D[rownames(Dtemp), ] <- as.matrix(Dtemp)
  dptemp <- dat.err$dosage.p[input.map$maps[[i.lpc]]$seq.num]
  dqtemp <- dat.err$dosage.q[input.map$maps[[i.lpc]]$seq.num]
  dq<-dp<-rep(m/2, length(mrknames))
  names(dp)<-names(dq)<-mrknames
  dp[names(dptemp)]<-dptemp
  dq[names(dqtemp)]<-dqtemp
  
  phP <- phQ <- vector("list", n.mrk)
  for(i in 1:length(phP)){
    phP[[i]] <- phQ[[i]] <- c(0:(m/2 - 1))
  }
  names(phP) <- names(phQ) <- mrknames
  phP[rownames(Dtemp)] <- input.map$maps[[i.lpc]]$seq.ph$P
  phQ[rownames(Dtemp)] <- input.map$maps[[i.lpc]]$seq.ph$Q
  
  ## Including error  
  gen<-vector("list", length(indnames))
  names(gen)<-indnames
  d.pq<-data.frame(dp = dp, 
                   dq = dq)
  d.pq$mrk<-rownames(d.pq)
  for(i in names(gen))
  {
    a<-matrix(0, nrow(D), input.map$info$m+1, dimnames = list(mrknames, 0:input.map$info$m))
    for(j in rownames(a)){
      if(D[j,i] == input.map$info$m+1){
        a[j,]<-segreg_poly(m = input.map$info$m, dP = dp[j], dQ = dq[j])
      } else {
        a[j,D[j,i]+1]<-1          
      }
    }
    a<-as.data.frame(a)
    a$mrk<-rownames(a)
    a.temp<-t(merge(a, d.pq, sort = FALSE)[,-c(1)])
    if(!is.null(error))
      a.temp<-apply(a.temp, 2, genotyping_global_error, 
                    m = input.map$info$m, error=error, 
                    th.prob = th.prob, 
                    restricted = restricted)
    else
      a.temp <- a.temp[1:(input.map$info$m+1), ]
    colnames(a.temp)<-a[,1]
    gen[[i]]<-a.temp
  }
  g = as.double(unlist(gen))
  p = as.numeric(unlist(phP))
  dp = as.numeric(cumsum(c(0, sapply(phP, function(x) sum(length(x))))))
  q = as.numeric(unlist(phQ))
  dq = as.numeric(cumsum(c(0, sapply(phQ, function(x) sum(length(x))))))
  rf = as.double(mf_h(diff(map.pseudo)))
  res.temp <-
    .Call(
      "calc_genoprob_prior",
      as.numeric(m),
      as.numeric(n.mrk),
      as.numeric(n.ind),
      as.numeric(p),
      as.numeric(dp),
      as.numeric(q),
      as.numeric(dq),
      as.double(g),
      as.double(rf),
      as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
      as.double(0),
      as.numeric(verbose),
      PACKAGE = "mappoly"
    )
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"),
                                mrknames, indnames)
  structure(list(probs = res.temp[[1]], map = map.pseudo), class="mappoly.genoprob")
}
meiotic_configuration<-function(a, pr.hom, parent){
  m = dim(pr.hom)[1]
  if(is.null(a[1])) return(NA)
  if(as.character(a[1])=="Zero") return(0)
  if(as.character(a[1])=="Inconclusive") return("i")
  # choosing parent
  if(parent == 1 ){
    lt<-letters[1:m]
  } else if(parent == 2){
    lt<-letters[(m+1):(2*m)]
  } else stop("parents should be '1' or '2'.")
  b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
  a$from <- sapply(b, function(x) x[1])
  a$to <- sapply(b, function(x) x[2])
  a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
  if(nrow(a) == 1) return(1) 
  x<-a[,3:4]
  g <- igraph::graph_from_data_frame(x)
  igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
  #plot(as.undirected(g))
  I<-igraph::components(g)$csize
  return(max(I))
#  if(all(I < 3)){
#    return("b")
#  } else if(any(I == 3) | any(I == 4)){
#    return("t")} else if(any(I > 4)) {
#      return("h")} else return("i") 
}
meiotic_graph<-function(a, pr.hom, parent){
  m = dim(pr.hom)[1]
  lt<-rownames(pr.hom)
  if(parent == 1){
    pal<-c("#800000", "#9A6324", "#808000", "#e6194B", "#f58231", "#ffe119", '#bcf60c', "#3cb44b")
    pal<-pal[ceiling(seq(1, 8, length.out = m))]
    names(pal)<-lt
  } else {
    pal<-c('#aaffc3',"#469990", "#42d4f4", "#4363d8", "#000075", "#911eb4", '#e6beff', '#f032e6')
    pal<-pal[ceiling(seq(1, 8, length.out = m))]
    names(pal)<-lt
  }
  if(is.null(a[1]))
  {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,m))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=rep(NA, m)))
  }
  if(as.character(a[1])=="Inconclusive")
  {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,m))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=rep(NA, m)))
  }  
  if(as.character(a[1])=="Zero") {
    g <- make_empty_graph() +
      vertices(lt)
    V(g)$size<-c(rep(40,m))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=pal))
  } else {
    b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
    a$from <- sapply(b, function(x) x[1])
    a$to <- sapply(b, function(x) x[2])
    a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
    x<-a[,3:4, drop = FALSE]
    remaining<-setdiff(lt, unique(as.character(as.matrix(x))))
    g <- igraph::graph_from_data_frame(d = x, directed = FALSE)
    if(length(remaining) > 0)
      g <- add_vertices(graph = g, nv = length(remaining), name=remaining)
    V(g)$size<-c(rep(40,m))
    E(g)$edge.color <- "gray80"
    E(g)$width <- 5
    V(g)$label.cex = 2
    igraph::V(g)$type <- substr(igraph::V(g)$name, 1, 1)=="V"
    return(list(g=g, vertex.color=pal[names(V(g))]))
  }
}
create_map_hap<- function(input.map, dat.dist, step = Inf,
                          phase.config = "best")
{
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  mrknames <- dat.dist$mrk.names[input.map$maps[[i.lpc]]$seq.num]
  map <- c(0, cumsum(imf_h(input.map$maps[[i.lpc]]$seq.rf)))
  names(map)<-mrknames
  if(is.null(step))
    return(map)
  minloc <- min(map)
  map <- map-minloc
  a <- seq(floor(min(map)), max(map), by = step)
  a <- a[is.na(match(a,map))]
  names(a) <- paste("loc",a,sep = "_")
  return(sort(c(a,map))+minloc)
}
multi_evidence<-function(res){
  if(as.character(res[1])=="Zero") return(FALSE)
  if(as.character(res[1])=="Inconclusive") return(FALSE)
  a<-res[!duplicated(res$`Hom. Pair.`),,drop = FALSE]
  b<-strsplit(as.character(a$`Hom. Pair.`), split = "")
  a$from <- sapply(b, function(x) x[1])
  a$to <- sapply(b, function(x) x[2])
  a[,c("from","to")]<-t(apply(a[,c("from","to")], 1, sort))
  return(any(duplicated(a$from)) | any(duplicated(a$to)))  
}

# individual <- "B2721.026"
# parent <- 2
# pr.hom <-  homo_prob(individual, genoprob, parent)
# hom.prob.thresh = 0.80
# seg.length.thresh = 10
# perc.info = 20
# dist.thresh = 2
# map.mappoly <- final.maps[[1]]
# dat.mappoly <- B2721
# map <- genoprob$map
# title.plot = NULL
# p <- count_cross(pr.hom, map, individual, parent, hom.prob.thresh, seg.length.thresh,
#                   perc.info, dist.thresh, map.mappoly, dat.mappoly, title.plot)
count_cross<-function(pr.hom, 
                      map,
                      individual,
                      parent,
                      hom.prob.thresh = 0.80,
                      seg.length.thresh = 10,
                      perc.info = 20,
                      dist.thresh = 2,
                      map.mappoly,
                      dat.mappoly,
                      title.plot = NULL)
{
  plot.1<-plot_homo_prob(pr.hom = pr.hom, map = map, ind = individual, title.plot = title.plot)
  pr.hom.2<-elim_low_prob_homo_segments(pr.hom = pr.hom, hom.prob.thresh = hom.prob.thresh)
  plot.2<-plot_homo_prob(pr.hom = pr.hom.2, map = map, ind = individual, title.plot = title.plot)
  segment.lengths<-get_segment_length(bool.pr.hom = pr.hom.2, map = map)
  pr.hom.3<-elim_small_and_incomplete_segments(bool.pr.hom = pr.hom.2, 
                                               segment.lengths = segment.lengths, 
                                               seg.length.thresh = seg.length.thresh)
  if(is.null(pr.hom.3)){
    pr.hom[]<-0
    return(list(p = "Inconclusive", 
                plot.1 = plot.1, 
                plot.final = NULL, 
                n.breaks=segment.lengths[[2]],
                n.breaks.new=NULL,
                res.par = "Inconclusive"))
  }
  if(100*ncol(pr.hom.3)/length(map) < perc.info){
    pr.hom[]<-0
    return(list(p = "Inconclusive", 
                plot.1 = plot.1, 
                plot.final = NULL, 
                n.breaks=segment.lengths[[2]],
                n.breaks.new=NULL,
                res.par = "Inconclusive"))
  }
  gap.lengths <- get_segment_length(bool.pr.hom = pr.hom.3, map = map, is.gap = TRUE)
  pr.hom.4<-elim_small_and_incomplete_segments(bool.pr.hom = pr.hom.3, 
                                               segment.lengths = gap.lengths, 
                                               seg.length.thresh = seg.length.thresh)
  if(is.null(pr.hom.4)){
    pr.hom[]<-0
    return(list(p = "Inconclusive", 
                plot.1 = plot.1, 
                plot.final = plot_homo_prob(pr.hom = pr.hom, map = map, ind = individual, title.plot = title.plot), 
                n.breaks=segment.lengths[[2]],
                n.breaks.new=NULL,
                res.par = "Inconclusive"))
  }
  if(100*ncol(pr.hom.4)/length(map) < perc.info){
    pr.hom[]<-0
    return(list(p = "Inconclusive", 
                plot.1 = plot.1, 
                plot.final = plot_homo_prob(pr.hom = pr.hom, map = map, ind = individual, title.plot = title.plot), 
                n.breaks=segment.lengths[[2]],
                n.breaks.new=NULL,
                res.par = "Inconclusive"))
  }
  ##plot.final <- plot_homo_prob(pr.hom = pr.hom.4, map = map, ind = ind, title.plot = title.plot)
  problematic.mrk<-setdiff(colnames(pr.hom), colnames(pr.hom.4))
  ## Reestimating missing segments using HMM
  if(100*(ncol(pr.hom) - length(problematic.mrk))/ncol(pr.hom) < perc.info){
    pr.hom[]<-0
    return(list(p = "Inconclusive",
                plot.1 = plot.1, 
                plot.final = plot_homo_prob(pr.hom = pr.hom, map = map, ind = individual, title.plot = title.plot), 
                n.breaks=segment.lengths[[2]], 
                res.par = "Inconclusive"))
  }
  dat.mappoly$geno.dose[problematic.mrk, individual] <- dim(pr.hom)[1] + 1  
  ## HMM estimation
  w<-calc_genoprob_err_one_ind(input.map = map.mappoly, dat.err = dat.mappoly, individual = individual, error = 0.05, verbose = FALSE)
  pr.hom.5 <- homo_prob(ind = individual, genoprob = w, parent = parent)
  p<-detect_co(pr.hom = pr.hom.5, map = map, dist.thresh = dist.thresh)
  plot.final <- plot_homo_prob(pr.hom = pr.hom.5, map = map, ind = individual, title.plot = title.plot)
  list(p = p$result,
       plot.1 = plot.1, 
       plot.final = plot.final, 
       n.breaks=segment.lengths[[2]], 
       res.par = p$result.par,
       pr.hom.final = pr.hom.5)
}

# individual <- "B2721.026"
# parent <- 2
# pr.hom <-  homo_prob(individual, genoprob, parent)
# hom.prob.thresh = 0.80
# seg.length.thresh = 10
# perc.info = 20
# dist.thresh = 2
# map.mappoly <- final.maps[[1]]
# dat.mappoly <- B2721
# map <- genoprob$map
# title.plot = NULL
# res <- plot_recombination_points(pr.hom, map, individual, parent, hom.prob.thresh = 0.8,
#                                  seg.length.thresh = 10, perc.info = 20,thresh.nbreak1 = Inf,
#                                  dist.thresh = 2, map.mappoly, dat.mappoly, title.plot = NULL)

plot_recombination_points<-function(pr.hom, map, individual, parent, hom.prob.thresh = 0.8, 
                                    seg.length.thresh = 10, perc.info = 20, thresh.nbreak1 = Inf, 
                                    dist.thresh = 2, map.mappoly, dat.mappoly, title.plot = NULL)
{
  p <- count_cross(pr.hom = pr.hom, 
                   map = map, 
                   individual = individual, 
                   parent = parent, 
                   hom.prob.thresh = hom.prob.thresh, 
                   seg.length.thresh = seg.length.thresh, 
                   dist.thresh = dist.thresh,
                   perc.info = perc.info,
                   map.mappoly = map.mappoly, 
                   dat.mappoly = dat.mappoly,
                   title.plot = title.plot)
  if(as.character(p[[1]][1])=="Inconclusive" | any(p$n.breaks > thresh.nbreak1)){
    return(list(plot=p[[3]], summary = p[[5]], pr.hom.final = p$pr.hom.final, 
                meiotic.configuration = "Inconclusive", 
                meiotic.graph = meiotic_graph(a = p[[1]], pr.hom = pr.hom, parent = parent)))
  } else if(as.character(p[[1]][1])=="Zero" | any(p$n.breaks > thresh.nbreak1)){
    return(list(plot=p[[3]], summary = p[[5]], pr.hom.final = p$pr.hom.final,  
                meiotic.configuration = meiotic_configuration(p[[5]], pr.hom = pr.hom, parent = parent), 
                meiotic.graph = meiotic_graph(a = p[[1]], pr.hom = pr.hom, parent = parent)))
  } else {
    p[[1]]<-p[[1]][nchar(as.character(p[[1]]$hom.pair))==2,]
    pb <- ggplot_build(p[[3]])
    pg <- ggplot_gtable(pb)
    for(i in unique(p[[1]]$hom.pair))
    {
      a<-subset(p[[1]], hom.pair%in%i)
      b<-pg$layout[2:7,]
      ystart <- 1/(2*(abs(a$from[1] - a$to[1]) + 1))
      yend <- 1 - ystart
      pnames1 <- names(pg$grobs[[a$from[1]+1]]$children)
      pname1 <- pnames1[str_detect(pnames1, "geom_density")]
      ####
      p1 <- pg$grobs[[a$from[1]+1]]$children[[pname1]]$children[[1]]$children[[1]]
      x1<-p1$x[1:(length(p1$x)/2)]
      pg <- gtable_add_grob(pg, segmentsGrob(x1[a$point.id], 
                                             rep(ystart, nrow(a)), 
                                             x1[a$point.id],
                                             rep(yend, nrow(a)), 
                                             gp = gpar(lty=1, lwd = 2, col = 1), 
                                             arrow = arrow(angle = 10, length = unit(0.15, "inches"), ends = "both", type = "closed")), 
                            t=b[a$to[1],"t"], b = b[a$from[1],"b"], l=b[a$from[1],"l"])
    }
    z1<-p[[5]]
    colnames(z1)<-c("Hom. Pair.", "CO position (cM)")
    return(list(plot=ggplotify::as.ggplot(pg), summary = z1, pr.hom.final = p$pr.hom.final, plot.orig = p$plot.1, 
                meiotic.configuration = meiotic_configuration(a = z1, pr.hom = pr.hom, parent = parent), 
                meiotic.graph = meiotic_graph(a = z1, pr.hom = pr.hom, parent = parent)))
  }
}

## For Potato paper
get_rec_chain <- function(x, parent){
  y <- numeric(5)
  names(y) <- c("0", "1", "2", "3", "4")
  if(parent == 1)
    z<-table(unlist(x$P1$c.p1))
  else
    z<-table(unlist(x$P2$c.p2))
  y[names(z)] <- z
  y
}



