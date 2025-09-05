# utils.R
# Internal utility functions for TExT

#' Identify sequential nodes in a phylogenetic tree
#'
#' Internal function to calculate which internal nodes in a phylogenetic tree 
#' have child lineages that branch sequentially (i.e., differ in rank by 1).
#'
#' @param tree An object of class \code{phylo}.
#'
#' @return An integer vector of length \code{2 * nTip - 1}, where the first 
#' \code{ntip} entries (corresponding to tips) are 0, and the remaining 
#' entries correspond to nodes. A value of 1 indicates a sequential node, 0 otherwise.
#'
#' @keywords internal
#'

snodes = function(tree) { 
  ntip = length(tree$tip.label)
  nd.edges =tree$edge[ which(tree$edge[,2] > ntip),]
  heights = node.depth.edgelength(tree)[(ntip+1):(2*ntip-1)]
  ranks = order(order(heights)) 
  rank.version.of.nd.edges = 0*nd.edges
  rank.version.of.nd.edges[,1]=ranks[nd.edges[,1]-ntip]
  rank.version.of.nd.edges[,2]=ranks[nd.edges[,2]-ntip]
  diffs = rank.version.of.nd.edges[,2]-rank.version.of.nd.edges[,1]
  node.ss = ifelse(diffs ==1, 1, 0)
  svals = c(0*(1:ntip), node.ss) 
  return(svals) 
} 

#' @param tree object of class phylo
#' @param m an integer; we ask whether each node's desc is among the next m 
#' @return A list. details: a data frame with success (was descn in the next m), 
#' probabilities (how probable was that), nodeids (as in tree); mW: W values for all 
#' offsets (eg 1, 3, 5, .. for m=2, and 2, 4, 6 ...). W: mean of mW. 
#' @keywords internal
descinm <- function(tree,m=2) {
  if (!inherits(tree, "phylo")) 
    stop("Tree should be an object of class \"phylo\".")
  lttb  <- lttbelow(tree)
  IND <- lttb$nodeids > length(tree$tip.label) 
  intonly=lttb$nodeids[IND] 
  nlins=lttb$lttbelow[IND]
  L=length(intonly)
  success=0*intonly
  probs=0*intonly 
  for (j in 1:(L-1)) {
    dj=tree$edge[which(tree$edge[,1]==intonly[j]),2]
    nextm=intonly[ (j+1):min(j+m,L)]
    nlj=nlins[j:min(j+m-1,L)] 
    success[j]=any(is.element(dj,nextm))
    probs[j]=prod(1 - 1/nlj)
  }
  mW=0*(1:m)
  for (Offset in 0:(m-1)) {
    Filter=seq(from=1+Offset,by=m,to=L-1) 
    mW[Offset+1]=(sum(success[Filter])) 
  }
  return(list(details=data.frame(s=success,probs=probs,nodeids=intonly),mW=mW,W=mean(mW)))
} 

################################################
#### are k of the next m to branch among the descendants of node i?
################################################
#' Are k of the next m to branch among the descendants of node i?
#' @param tree object of class phylo
#' @param k integer: are k of next m descs descended from node i?
#' @param m integer: are k of next m descs descended from node i?
#' @return A list.  details: a data frame with success (k desc of i in the next m), 
#' probabilities (how probable was that), nodeids (as in tree); kmW: W values for all 
#' offsets (eg 1, 3, 5, .. for m=2, and 2, 4, 6 ...). W: mean of kmW.
#' @keywords internal
kinm <- function(tree,k=2,m=3) {
  if (!inherits(tree, "phylo")) 
    stop("Tree should be an object of class \"phylo\".")
  num.tips=length(tree$tip.label)
  lttb<-lttbelow(tree)
  IND<-lttb$nodeids > length(tree$tip.label) 
  intonly=lttb$nodeids[IND] 
  nlins=lttb$lttbelow[IND] 
  L=length(intonly) 
  mck=choose(m,k)
  success=0*intonly 
  probs=0*intonly
  for ( j in 1:(L-m)) {
    dj=phytools::getDescendants(tree,intonly[j])
    Nj=(sum(dj <= num.tips) -2) 
    nextm=intonly[ (j+1):min(j+m,L)] 
    success[j]=sum(is.element(dj,nextm))>=k
    probs[j]= choose(Nj,k)*choose( L-j-Nj, m-k)/choose(L-j,m)
  }
  mW=0*(1:m)
  for (Offset in 0:(m-1)) {
    Filter=seq(from=1+Offset,by=m,to=L-1)
    mW[Offset+1]=(sum(success[Filter]))
  } 
  return(list(details=data.frame(success=success,probs=probs,nodeids=intonly),kmW=mW,W=mean(mW)))
}

#' Compute times of nodes and tips and arrange in chronological order along with LTT
#' @param tree Object of class phylo, or convertible with as(tree,"phylo")
#' @return data frame with times, topdepths (topological depth from root, in number of edges, 
#' lttbelow: the number of lineages in the tree just below the current node, and nodeids: the ID of the node. 
#' 1: number tips are the tips. Order corresponds to the order in the phylo tree. 
#' @keywords internal 
lttbelow <- function(tree) {
  if (!inherits(tree,"phylo")) {tree=as(tree,"phylo")}
  N=2*tree$Nnode+1 
  Ntips=(N+1)/2 
  times=vector(mode="numeric",length=N)
  times[Ntips+1]=0; 
  topdepths=times 
  for (k in 1:(N-1)) {
    times[tree$edge[k,2]]=tree$edge.length[k]+ times[tree$edge[k,1]]
    topdepths[tree$edge[k,2]]=1+topdepths[tree$edge[k,1]]
  }
  nodeids=order(times)
  sorttimes=times[nodeids]
  contribs=-1+2*as.numeric(order(times)>Ntips) # each contribution to ltt
  ltt=1+cumsum(contribs)
  return(data.frame(times=sorttimes,topdepths=topdepths[nodeids], lttbelow=ltt,nodeids=nodeids))
}

#' Compute times of nodes and tips and arrange in chronological order 
#' along with LTT, if the tree has polytomies. 
#' @param tree Object of class phylo, or convertible with as(tree,"phylo")
#' @return data frame with times, topdepths (topological depth from root, in number of edges, 
#' lttbelow: the number of lineages in the tree just below the current node, and nodeids: the ID of the node. 
#' 1: number tips are the tips. Order corresponds to the order in the phylo tree.
#' @keywords internal
lttpoly <- function(tree) {
  if (!inherits(tree,"phylo")) tree=as(tree,"phylo")
  N = tree$Nnode + length(tree$tip.label) 
  Ntips=length(tree$tip.label)
  times=vector(mode="numeric",length=N)
  times[Ntips+1]=0; 
  topdepths=times 
  for (k in 1:(N-1)) {
    times[tree$edge[k,2]]=tree$edge.length[k]+ times[tree$edge[k,1]]
    topdepths[tree$edge[k,2]]=1+topdepths[tree$edge[k,1]]
  }
  nodeids=order(times) 
  sorttimes=times[nodeids]
  nums.desc = vapply(nodeids, function(x) getnumdesc(tree,x), FUN.VALUE = 1)
  contribs=-1+nums.desc 
  ltt=1+cumsum(contribs)
  return(data.frame(times=sorttimes,topdepths=topdepths[nodeids], 
                    lttbelow=ltt,nodeids=nodeids))
}

#' Get the number of child lineages descending from a node
#' @param tree An object of class \code{phylo}.
#' @param node Integer; the node index in \code{tree}.
#' @return An integer giving the number of direct child lineages of the node. 
#' For tips, returns 0.
#' @keywords internal
getnumdesc = function(tree, node) { 
  if (node <= length(tree$tip.label)) {
    return(0)
  } else { 
    return(length(which(tree$edge[,1]==node))) 
  }
} 

#' @description find the date of the root of a phylogeny
#' @param tree a phylogenetic tree in phylo format 
#' @param info a data frame with dates ('Date') and tip labels (TipLabels_Data)
#' @return the date of the root
#' @keywords internal
getrootdate <- function(tree, info) {
  myind = sample(1:length(tree$tip.label), 1) 
  mytip = tree$tip.label[myind] 
  tdate =  info$Date[which(info$TipLabels_Data == mytip)]
  tdepth = node.depth.edgelength(tree)[myind]
  rootdate = tdate - days(floor(365*tdepth))
  return(rootdate)
}

#' @param tree Object of class phylo
#' @return Data frame linking tip names to sequential branching along paths to root
#' @keywords internal
tipprofiles <- function(tree) {
  if (!inherits(tree,"phylo")) {tree=as(tree,"phylo")}
  N=length(tree$tip.label)
  sstuff=getstattest(tree)$details
  sinorder=sstuff$s[order(sstuff$nodeids)]
  lttstuff=lttbelow(tree)
  rightorder=order(lttstuff$nodeids)
  timesinorder=lttstuff$times[rightorder[-(1:N)]]
  depthsinorder=lttstuff$topdepths[rightorder[-(1:N)]]
  timeprod=timesinorder*sinorder 
  topdepprod=depthsinorder*sinorder
  tipsumofs=vector(mode="numeric",length=N)
  tipsumtimeprod=tipsumofs; tipsumtdprod=tipsumofs
  sumofs=vector(mode="numeric",length=N-1)
  stimeprod=sumofs; stopdepprod=sumofs 
  for (k in 1:(N-1)) {
    if (tree$edge[k,2]>N) {
      sumofs[tree$edge[k,2]-N]=sumofs[tree$edge[k,1]-N]+sinorder[tree$edge[k,2]-N]
      stimeprod[tree$edge[k,2]-N]=stimeprod[tree$edge[k,1]-N]+timeprod[tree$edge[k,2]-N]
      stopdepprod[tree$edge[k,2]-N]=stopdepprod[tree$edge[k,1]-N]+topdepprod[tree$edge[k,2]-N]
    } else {
      tipsumofs[tree$edge[k,2]]=sumofs[tree$edge[k,1]-N]
      tipsumtimeprod[tree$edge[k,2]]=stimeprod[tree$edge[k,1]-N]
      tipsumtdprod[tree$edge[k,2]]=stopdepprod[tree$edge[k,1]-N]                
    }
  }
  FIXIND <- order(tree$tip.label)
  return(data.frame(SumS=tipsumofs[FIXIND], sWithLengths=tipsumtimeprod[FIXIND],sWithTopDepths=tipsumtdprod[FIXIND],TipNames=tree$tip.label[FIXIND]))
}