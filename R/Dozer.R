#' Compute gene-gene correlation matrix
#'
#' This function take gene expression matrix (gene by cell), cell sequencing depth 
#' and other cell level covariates that need to be accounted for in expression normalization as input to 
#' compute gene-gene correlation matrix.
#' @param data A gene expression matrix, with rows representing genes and columns representing cells.
#' @param lib_size A vector of cell sequencing depth (optional).
#' @param covs A data frame or a list of other covariates for cells, e.g. batch labels.
#' @return co-expression matrix in the "network" slot and gene noise to signal ratio vector in the "ratio" slot.
#' @export
compute_gene_correlation <- function(data, lib_size = NULL, covs = NULL){
  ## shape of data
  nr = nrow(data)
  nc = ncol(data)
  ## If cell sequencing depth is not provided, substitute it with total counts per cell
  if (is.null(lib_size)){
    lib_size = colSums(data)
  }
  ## Re-scale lib_size in case it is too large of small in magnitude
  lib_size = lib_size/median(lib_size)
  
  dat_normalize = .normalize_data(data = data, lib_size = lib_size , covs = covs)
  cell_weight = lib_size
  
  ## Re-scale the normalized counts to avoid numerical instability when computing correlations
  norm_data = t(apply(dat_normalize$norm_expr, 1, FUN=function(x){x/mean(x)}))
  raw_corr = cov.wt(t(norm_data), wt = cell_weight, cor = TRUE)$cor
  diag(raw_corr) = 0
  
  scaling = .noise_ratio_correction(dat_normalize$norm_expr, dat_normalize$expr_var, cell_weight)
  network = sweep(raw_corr, MARGIN=2, scaling$correction_factor, `*`)
  network = sweep(network, MARGIN=1, scaling$correction_factor, `*`)
  rownames(network) = rownames(data)
  colnames(network) = rownames(data)
  ratio = data.frame(ratio = scaling$noise_ratio)
  rownames(ratio) = rownames(data)
  return(list(network = network, ratio = ratio))
}

#' Compute gene centrality in the hard-threshold gene co-expression network.
#'
#' @param network Gene-gene correlation matrix.
#' @param threshold quantile threshold on the absolute correlations for hard-thresholding.
#' @return list of gene centrality vectors, including degree, pagerank, betweenness and eigenvector centrality.
#' @export
compute_centrality <- function(network, threshold = 0.99){
  diag(network) = 0
  network = abs(network)
  network[network < quantile(network[upper.tri(network)], threshold)]=0
  network[network > 0] = 1
  node = nrow(network)
  pagerank = rep(0, node)
  degree = rep(0, node)
  betweenness = rep(0, node)
  ev = rep(0, node)
  
  igraph_net = igraph::graph_from_adjacency_matrix(network, mode =  "undirected")
  pagerank =  igraph::page_rank(igraph_net)$vector
  betweenness = igraph::betweenness(igraph_net, directed = F)
  degree = igraph::degree(igraph_net)
  ev = igraph::evcent(igraph_net)$vector
  
  return(list(degree = degree, 
              pagerank = pagerank, 
              betweenness = betweenness, 
              eigenvector = ev))
}

#' Simulate gene expression matrix from a Gamma-Poisson distribution.
#'
#' @param corr Gene-gene correlation matrix served as ground truth.
#' @param shape,scale Gamma shape and scale parameter in the Gamma-Poisson distribution.
#' @param mlog,sdlog Mean and standard deviation of the log normal distribution for library size.
#' @param ncell Number of cells to be simulated.
#' @return simulated gene expression matrix using Gamma Poisson distribution.
#' @export
simulate_counts<- function(corr, shape, scale, mlog, sdlog, ncell){
  ngene = nrow(corr)
  G = t(lcmix::rmvgamma(ncell, shape = shape, rate = 1, corr = corr))
  lib0 = rlnorm(ncell, meanlog = mlog, sdlog = sdlog) 
  nij=(matrix(scale, nrow=ngene, ncol=1)%*% matrix(lib0, nrow = 1))
  X=G*nij
  data = matrix(rpois(ngene* ncell, lambda = as.vector(X)), nrow=ngene)
  rownames(data)= paste0('gene', 1:ngene)
  colnames(data)=paste0('cell', 1:ncell)
  return(data)
}

#' Conduct hierarchical clustering for "difference network" between diagnosis group.
#'
#' @param network1 network in diagnosis group 1
#' @param network2 network in diagnosis group 2
#' @return cluster labels
#' @export
clustering_difference_network<-function(network1, network2, minClusterSize = 20){
  diff_network = network1 - network2
  # compute inner product of difference network to itself, 
  # so that genes with similar connectivity patterns in the difference network
  # have high score in the transformed adjacency matrix
  diff_network = diff_network%*%t(diff_network)
  
  non_singleton = rowMeans(abs(diff_network))>=0
  sub_net = diff_network[non_singleton, non_singleton]
  ## convert to a nonnegative distance matrix
  distance = max(sub_net) - sub_net
  cluster_non_singleton = dynamicTreeCut::cutreeDynamic(dendro = stats::hclust(stats::as.dist(distance)), 
                                      distM = distance, minClusterSize = minClusterSize)
  final_cluster = rep('singleton', nrow(network1))
  final_cluster[non_singleton] = cluster_non_singleton
  names(final_cluster) = rownames(network1)
  return(final_cluster)
}



## internal functions ##
.normalize_data <- function(data, lib_size, covs = NULL){
  # a function that uses poisson regression to normalize count matrix of scRNA-seq data 
  # input: data (count matrix), lib_size (a vector of cell sequencing depth), covs (other sources of variation to be adjusted)
  # output a list containing norm_expr (normalized matrix) and expr_var (variance of the normalized count)
  ngene = nrow(data)
  ncell = ncol(data)
  norm_expr = matrix(0, nrow = ngene, ncol = ncell)
  expr_var = matrix(0, nrow = ngene, ncol = ncell)
  
  cores = parallel::detectCores()
  cl <- parallel::makeCluster(cores[1]-1) 
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`
  res = foreach::foreach(i = 1:ngene, .combine = 'rbind')%dopar%{
    if (length(covs)==0){
      obs = data.frame(y = data[i,], n = lib_size)
    }else{
      obs = data.frame(y = data[i,], n = lib_size, covs)
    }
    
    fit1 = glm(y  ~ .-n, offset=(log(n)),
               family=poisson(link = log), data = obs)
    res_var = obs$y/(obs$n)^2
    res = resid(fit1, type = "response")/obs$n
    c(res, res_var)
  }
  parallel::stopCluster(cl)
  return(list(norm_expr = res[,1:ncell], expr_var = res[,(ncell+1):(ncell*2)]))
  
}


.noise_ratio_correction <- function(norm_expr, expr_var, cell_weight){
  # a function that computes correction factor for each gene
  # input: norm_expr (normalized matrix), expr_var (variance of the normalized count), cell_weight
  # output: gene noise ratio and correction factor
  ncell = ncol(norm_expr)
  l = apply(norm_expr, 1, FUN=function(x){matrixStats::weightedVar(x, cell_weight)})
  e = apply(expr_var, 1, FUN=function(x){weighted.mean(x, cell_weight)}) 
  x1 = t(apply(norm_expr, 1, FUN=function(x){cell_weight*(x - weighted.mean(x, cell_weight))^2}))
  x2 = t(apply(expr_var, 1, FUN=function(x){cell_weight*x}))
  a = matrixStats::rowVars(x1)/ncell
  b = matrixStats::rowVars(x2)/ncell
  
  v = pmax(0,(a*e^2+b*l^2)/(l-e)^4)
  y = pmax(1,l/(l-e)/(1+v))    
  x = e/l
  z = loess(y~x, span = 0.1, degree = 1)
  m = max(z$fitted)
  y[x >= 1-1/m] = m
  y[x<1-1/m] = 1/(1-x[x<1-1/m])
  
  return(list(noise_ratio = e/l, correction_factor = sqrt(y)))
}


