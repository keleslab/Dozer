#' Compute gene-gene correlation matrix
#'
#' This function take gene expression matrix (gene by cell), cell sequencing depth 
#' and other cell level covariates that need to be accounted for in expression normalization as input to 
#' compute gene-gene correlation matrix.
#' @param data A gene expression matrix, with rows representing genes and columns representing cells.
#' @param lib_size A vector of cell sequencing depth (optional).
#' @param covs A data frame, whose columns are covariates for cells, e.g. batch labels.
#' @param multicore A boolean variable indicating whether we want to parallel the computation for the normalization of each gene.
#' @param min_expressed_cells The minimum number of cells that a gene have positive expression in. The co-expression and noise ratio of genes expressed in less than this number of cells will not be computed and set as NA (missing).
#' @param gene_group_quantile A vector of quantile values (ranging from 0 to 1) on gene average expression, to divide genes into groups for the estimation of cell size specifically of each gene.
#' @return co-expression matrix in the "network" slot and gene noise to signal ratio vector in the "ratio" slot.
#' @export
compute_gene_correlation <- function(data, lib_size = NULL, covs = NULL, multicore = FALSE, min_expressed_cells = 2, 
                                     gene_group_quantile = NULL){
  # convert sparse matrix to dense matrix
  if (inherits(data,'Matrix')){
    data = matrix(data, nrow = nrow(data), dimnames = list(rownames(data), colnames(data)))
  }
  ## store gene names
  if (is.null(rownames(data))){
    rownames(data) = paste0('Gene_', 1:nrow(data))
  }
  
  ## Give higher cells weights to cells of higher sequencing depth, in the computation of gene noise ratio and weighted gene correlation.
  cell_weight = colSums(data)/median(colSums(data))  
  
  all_genes = rownames(data)
  ## detect sparse genes
  sparse_gene = rowSums(data > 0) < min_expressed_cells
  genes = all_genes[!sparse_gene]
  
  ## If cell sequencing depth is not provided, substitute it with total counts per cell
  if (is.null(lib_size)){
    lib_size = trimmed_total_counts(data)
  }
  ## Re-scale lib_size to avoid numerical instability
  lib_size = lib_size/median(lib_size)
  
  if (!is.null(gene_group_quantile)){
    if (is.null(covs)){
      covs = .compute_tUMI_for_gene_bins(data[!sparse_gene, ], gene_group_quantile)
    }else{
      covs = cbind(covs, .compute_tUMI_for_gene_bins(data[!sparse_gene, ], gene_group_quantile))
    }
    lib_size = rep(1, ncol(data))
  }
  ## Normalize raw counts and compute gene noise ratio using internal function ".normalize_data"
  dat_normalize = .normalize_data(data = data[!sparse_gene,], lib_size = lib_size , covs = covs, multicore = multicore)
  
  ## Compute gene noise ration and gene correction factor using internal function ".noise_ratio_correction"
  scaling = .noise_ratio_correction(dat_normalize$norm_expr, dat_normalize$expr_var, cell_weight)
  
  ## Re-scale the normalized counts to avoid numerical instability when computing correlations
  dat_normalize$norm_expr = t(apply(dat_normalize$norm_expr, 1, FUN=function(x){x/mean(x)}))
  ## Compute weighted gene correlations.
  network_nonsparse_gene = cov.wt(t(dat_normalize$norm_expr), wt = cell_weight, cor = TRUE)$cor
  
  ## Apply gene correction factors 
  network_nonsparse_gene = sweep(network_nonsparse_gene, MARGIN = 2, scaling$correction_factor, `*`)
  network_nonsparse_gene = sweep(network_nonsparse_gene, MARGIN = 1, scaling$correction_factor, `*`)
  ## Label each row and column of gene correlation matrix with gene id.
  rownames(network_nonsparse_gene) = genes
  colnames(network_nonsparse_gene) = genes
  
  ## Set the correlation of a gene to itself to be 1
  diag(network_nonsparse_gene) = 1
  ## Threshold all corrected correlation larger than 1 to be the largest correlation value which is smaller than 1
  network_nonsparse_gene[network_nonsparse_gene > 1] = max(network_nonsparse_gene[network_nonsparse_gene < 1])
  
  if (sum(sparse_gene)>0){
    ## If there are sparse genes, set the correlation and noise ratio of these genes as missing (NA).
    network = matrix(NA, nrow = length(all_genes), ncol = length(all_genes))
    rownames(network) = all_genes
    colnames(network) = all_genes
    network[genes, genes] = network_nonsparse_gene
    ratio = data.frame(matrix(NA, nrow = length(all_genes), ncol = 2))
    colnames(ratio) = c('noise_ratio', 'correction_factor')
    rownames(ratio) = all_genes
    ratio[genes, 1] = scaling$noise_ratio
    ratio[genes, 2] = scaling$correction_factor
    
    return(list(network = network, ratio = ratio))
  }else{
    ## Label gene noise ratio with gene id.
    ratio = data.frame(scaling)
    rownames(ratio) = all_genes
    ## Return gene co-expression matrix and gene noise ratio.
    return(list(network = network_nonsparse_gene, ratio = ratio))
  }
  
  
}

#' Compute gene centrality in the hard-threshold gene co-expression network.
#'
#' @param network Gene-gene correlation matrix.
#' @param threshold quantile threshold on the absolute correlations for hard-thresholding.
#' @return list of gene centrality vectors, including degree, pagerank, betweenness and eigenvector centrality.
#' @export
compute_centrality <- function(network, threshold = 0.99){
  ## Set the diagonal of gene correlation matrix as 0, because we do not want network edges of genes to themselves.
  diag(network) = 0
  ## Take the absolute value of correlations to build an un-signed network.
  network = abs(network)
  ## Correlation on the quantile threshold.
  quant_value = quantile(network[upper.tri(network)], threshold)
  ## Set the gene pairs with correlation smaller than quant_value as 0 and 1 otherwise.
  network[network < quant_value] = 0
  network[network > 0] = 1
  
  ## Transform adjacency matrix into igraph (dependent package) format and then compute gene centrality using igraph functions.
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
simulate_counts <- function(corr, shape, scale, mlog, sdlog, ncell){
  ## The number of genes matches the dimension of correlation matrix "corr".
  ngene = nrow(corr)
  ## Simulate a matrix of ngene * ncell from multivariate Gamma distribution with shape = shape and rate = 1, and correlation corr.
  G = t(lcmix::rmvgamma(ncell, shape = shape, rate = 1, corr = corr))
  ## Multiply scale parameter of each row of matrix G, by the scaling property of Gamma, the Gamma scale parameter will equal to variable scale. 
  G = sweep(G, MARGIN = 1, scale, `*`)
  ## Simulate library size from a log normal distribution with parameter mlog and sdlog and multiply it to columns of G.
  lib0 = rlnorm(ncell, meanlog = mlog, sdlog = sdlog) 
  G = sweep(G, MARGIN = 2, lib0, `*`)
  ## Simulate the final Gamma-Poisson counts using the multivariate gamma G as Poisson rate.
  data = matrix(rpois(ngene* ncell, lambda = as.vector(G)), nrow = ngene)
  ## Give names to the row and columns of simulated counts and return.
  rownames(data) = paste0('gene', 1:ngene)
  colnames(data) = paste0('cell', 1:ncell)
  return(data)
}

#' Conduct hierarchical clustering for "difference network" between diagnosis group.
#'
#' @param network1 network in diagnosis group 1
#' @param network2 network in diagnosis group 2
#' @return cluster labels
#' @export
clustering_difference_network<-function(network1, network2, minClusterSize = 20){
  if (is.null(rownames(network1))){
    genes = paste0('Gene_', 1:nrow(network1))
  }else{
    genes = rownames(network1)
  }
  ## Compute the difference of the two adjacency matrices.
  diff_network = network1 - network2
  ## compute inner product of difference network to itself, 
  ## so that genes with similar connectivity patterns in the difference network
  ## have high score in the transformed adjacency matrix
  diff_network = diff_network%*%t(diff_network)
  
  ## separate the genes with edges in diff_network (non_singleton) and genes without edges in diff_network (singleton).
  non_singleton = rowSums(abs(diff_network)) > 0
  sub_net = diff_network[non_singleton, non_singleton]
  ## Convert the adjacency matrix of all non_singletons to a nonnegative distance matrix.
  distance = max(sub_net) - sub_net
  ## Conduct hierarchical clustering and dynamic Tree Cut on the non_singleton sub-network.
  cluster_non_singleton = dynamicTreeCut::cutreeDynamic(dendro = stats::hclust(stats::as.dist(distance)), 
                                                        distM = distance, minClusterSize = minClusterSize)
  ## Combine cluster labels for singleton and nonsingletons, labels genes by their id provided in network1.
  final_cluster = rep('singleton', nrow(network1))
  final_cluster[non_singleton] = cluster_non_singleton
  names(final_cluster) = genes
  return(final_cluster)
}

#' Compute trimmed total counts
#'
#' @param data counts matrix
#' @param max_weight genes with weights higher than max_weight will be trimmed 
#' @param delta A threshold on the minimum total counts over median total counts, in order not to have extremely small trimmed total counts, which aims to stablize normalized counts.
#' @return trimmed total counts per cell
#' @export
trimmed_total_counts <- function(data, max_weight = 1/300, delta = NULL){
  if (inherits(data,'Matrix')){
    data = matrix(data, nrow = nrow(data))
  }
  if(is.null(delta)){
    # raw total counts
    raw = colSums(data)
    # set a default value for delta
    delta = max(min(raw)/median(raw)/10, 1e-2)  
  }
  # Do not use extra sparse genes (less than 10 positive counts) in the calculation of trimmed UMI counts.
  non_sparse = rowSums(data > 0) > 10
  if (sum(non_sparse)<50){warning('Majority of genes are overly sparse or not enough genes provided in the dataset.')}
  data = data[non_sparse, ]
  
  N = nrow(data)
  max_weight = max(max_weight, 1/N)
  
  # Find the cutoff value for trimmed UMI through bisection
  weight = rowSums(data)
  weight = weight/sum(weight)
  if (max_weight==1/N){
    weight = rep(1/N, N)
  }else if(max(weight)> max_weight){
    s0=sort(weight)
    l = 1
    r = N
    while (l<r){
      s = s0
      m = floor(l/2+r/2)
      s[m:N]=s[m]
      M = max(s)/sum(s)
      if (M<max_weight){
        l = m + 1
      }else if(M==max_weight){
        l = m
        break
      }else{
        r = m - 1
      }
    } 
    weight[weight>=s0[l]] = s0[l]
    weight = weight/sum(weight)
  }
  
  # compute total UMI counts with gene weights specified by "weight"
  ans = rep(0, ncol(data))
  for (i in 1:nrow(data)){
    ans = ans+weight[i]*data[i,]/mean(data[i,])
  }
  ans = ans/median(ans)
  ans[ans<delta] = delta
  return(ans)
}

#' Generate a figure showing the slope of gene counts over trimmed total UMI counts.
#' If all gene slopes are centered around 1, using trimmed total UMI counts as adjustment for cell size is adequate.
#' @param data counts matrix
#' @param n number of bins for genes 
#' @return density plot of linear regression slopes of gene expression over trimmed total UMI counts.
#' @export
diagnoistic_plot_cell_size <- function(data, normalized_data = NULL, gene_group_quantile = NULL, n = 10){
  if (inherits(data,'Matrix')){
    data = matrix(data, nrow = nrow(data))
  }
  mean_expr = rowMeans(data)
  data_scale_row = t(apply(data, 1, FUN=function(x){x/mean(x)}))
  tUMI = trimmed_total_counts(data)
  tUMI = tUMI/mean(tUMI)
  df_raw = data.frame(expr = mean_expr, 
                      slope = apply(data_scale_row, 1, FUN=function(x){lm(x~tUMI-1)$coef[1]}), 
                      data = 'Raw counts')
  df_raw$expr_bins = ggplot2::cut_number(df_raw$expr, n=n)
  p1 = ggplot2::ggplot(df_raw, aes(slope, color= expr_bins))+ggplot2::scale_color_viridis_d()+
    ggplot2::geom_density(size=2) + ggplot2::geom_vline(xintercept = 1) + theme_bw(base_size = 10)+
    xlab('Slope of Raw counts\non trimmed total UMI counts')+labs(color = 'Gene bins\nby expression')
  
  if (is.null(normalized_data)){
    if (is.null(gene_group_quantile)){
      normalized_data = .normalize_data(data, tUMI)$norm_expr
    }else{
      covs = .compute_tUMI_for_gene_bins(data, gene_group_quantile)
      normalized_data = .normalize_data(data, lib_size = rep(1, ncol(data)), covs = covs)$norm_expr
    }
  }
  
  #normalized_data = t(apply(normalized_data, 1, FUN=function(x){x/mean(x)}))
  df_norm = data.frame(expr = rowMeans(data), 
                       cor = apply(normalized_data, 1, FUN=function(x){cor(x,tUMI)}), 
                       data = 'Normalized data')
  df_norm$expr_bins = ggplot2::cut_number(df_norm$expr, n=n)
  p2 = ggplot2::ggplot(df_norm, aes(cor, color= expr_bins))+ggplot2::scale_color_viridis_d()+
    ggplot2::geom_density(size=2) + ggplot2::geom_vline(xintercept = 0) + theme_bw(base_size = 10)+
    xlab('Correlation of normalized genes\nwith trimmed total UMI counts')+labs(color =  'Gene bins\nby expression')+
    xlim(-1,1)
  
  return(list(raw_data = p1, normalized_data = p2))
}



## internal functions ##
.compute_tUMI_for_gene_bins <- function(data, gene_group_quantile){
  # generate a data frame, whose columns are log of trimmed total UMI counts for genes in bins set by expression quantiles from "gene_group_quantile".
  nbins = length(gene_group_quantile)+1
  covs = data.frame(matrix(0, nrow = ncol(data), ncol = nbins))
  mean_expr = rowMeans(data)
  q = c(0, unlist(lapply(gene_group_quantile, FUN=function(x){quantile(mean_expr, x)})), max(mean_expr))
  for(i in 1:nbins){
    covs[,i] = log(trimmed_total_counts(data[mean_expr > q[i]& mean_expr <= q[i+1],]))
  }
  return(covs)
}

.normalize_data <- function(data, lib_size, covs = NULL, multicore = FALSE, delta = NULL){
  # a function that uses poisson regression to normalize count matrix of scRNA-seq data 
  # input: data (count matrix, gene by cell), lib_size (a vector of cell sequencing depth), covs (other sources of variation to be adjusted)
  # output a list containing norm_expr (normalized matrix) and expr_var (variance of the normalized count)
  
  ## shape of data is # genes by # cells.
  ncell = ncol(data)
  ngene = nrow(data)
  if (is.null(delta)){
    total_counts = colSums(data)
    delta = max(min(total_counts)/median(total_counts), 0.01)
  }
  
  
  ## remove covariates with identical values for all cells
  if (!is.null(covs)){
    covs = covs[, apply(covs, 2, FUN=function(i){length(unique(i))})!=1]
  }
  ## The number of covariate is zero, scale raw counts by cell lib_size, otherwise run a Poisson regression to determine cell specific size factor
  if (length(covs)==0){
    ## When the number of covariate is zero, scale raw counts by cell lib_size, 
    ## the variance of the normalized count is count/lib_size^2.
    lib_size = lib_size/median(lib_size)
    lib_size[lib_size<delta] = delta
    norm_expr = sweep(data, MARGIN = 2, 1/lib_size, `*`)
    expr_var = sweep(data, MARGIN = 2, 1/lib_size^2, `*`)
    return(list(norm_expr = norm_expr, expr_var = expr_var))
  }else if (!multicore){
    ## If there are covariates and we decide not to parallelize the computation for each gene,
    ## we run a Poisson regression for each gene in a loop and compute the response prediction as cell size factor for normalization.
    norm_expr = matrix(0, nrow = ngene, ncol = ncell)
    expr_var = matrix(0, nrow = ngene, ncol = ncell)
    for(i in 1:ngene){
      obs = data.frame(y = data[i,], n = lib_size, covs)
      fit1 = stats::glm(y  ~ .-n, offset=(log(n)),
                        family=poisson(link = log), data = obs)
      lhat = stats::predict(fit1, type='response')
      lhat = lhat/median(lhat)
      lhat[lhat<delta] = delta
      expr_var[i,] = obs$y/lhat^2
      norm_expr[i,] = obs$y/lhat
    } 
    return(list(norm_expr = norm_expr, expr_var = expr_var))
  }else{
    ## The same computation as the previous block with parallelization for each gene.
    cores = parallel::detectCores()
    cl <- parallel::makeCluster(cores[1]-1) 
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    res = foreach::foreach(i = 1:ngene, .combine = 'rbind')%dopar%{
      obs = data.frame(y = data[i,], n = lib_size, covs)
      fit1 = stats::glm(y  ~ .-n, offset=(log(n)),
                        family=poisson(link = log), data = obs)
      lhat = stats::predict(fit1, type='response')
      lhat = lhat/median(lhat)
      lhat[lhat<delta] = delta
      
      norm_expr_i = obs$y/lhat
      expr_var_i = obs$y/lhat^2
      c(norm_expr_i, expr_var_i)
    }
    parallel::stopCluster(cl)
    return(list(norm_expr = res[,1:ncell], expr_var = res[,(ncell+1):(ncell*2)]))
  }
}



.noise_ratio_correction <- function(norm_expr, expr_var, cell_weight){
  # A function that computes noise ration and correction factor for each gene.
  # Input: norm_expr (normalized matrix), expr_var (variance of the normalized count), cell_weight.
  # Output: gene noise ratio and correction factor.
  
  ## Shape of norm_expr is number of genes by number of cells.
  ncell = ncol(norm_expr)
  ngene = nrow(norm_expr)
  
  ## Dividend and divisor of supplement equation (7).
  mu_w = apply(expr_var, 1, FUN=function(x){weighted.mean(x, cell_weight)}) 
  s2_w = apply(norm_expr, 1, FUN=function(x){matrixStats::weightedVar(x, cell_weight)})
  ## Noise ratio.
  noise_ratio = mu_w/s2_w
  ## Because the dividend and divisor are  all estimators (with noise), the ratio can be slightly larger than 1.
  ## If noise ratio is larger than 1, the plug-in estimator for gene correction will be negative (supplement equation (10)).
  ## So we remove the genes with noise ratio greater than 1 in the smoothing procedure.
  noise_ratio_less_1 = (mu_w < s2_w)
  
  ## Estimate variance of the estimator \sqrt{\hat{S}_j}, supplement equation (11)
  s2_w2 = t(apply(norm_expr, 1, FUN=function(x){cell_weight*(x - weighted.mean(x, cell_weight))^2}))
  mu_w2 = t(apply(expr_var, 1, FUN=function(x){cell_weight*x}))
  var_s2_w = matrixStats::rowVars(s2_w2)/ncell
  var_mu_w = matrixStats::rowVars(mu_w2)/ncell
  var_s = (var_s2_w*mu_w^2+var_mu_w*s2_w^2)/(s2_w-mu_w)^4
  #var_root_s = (s2_w*var_mu_w + var_s2_w*mu_w^2/s2_w)/4/(s2_w-mu_w)^3 
  
  ## Shrinkage estimator S_j^0 = 1 v (\hat{S}_j)/(1+var(\hat{S}_j)) supplement equation (12) 
  S_shrinkage = pmax(1, 1/(1-noise_ratio[noise_ratio_less_1])/(1 + var_s[noise_ratio_less_1]))    
  
  ## LOESS smoothing for s_shrinkage
  S_smooth = loess(S_shrinkage ~ noise_ratio[noise_ratio_less_1], span = 0.1, degree = 1)
  ## Find the turning point of the uni-modal smoothing curve and truncate the correction factor at the turning point. (supplement equation (13))
  turn_val = max(S_smooth$fitted)
  turn_arg =  S_smooth$x[which.max(S_smooth$fitted)]#1-1/turn_val
  
  final_S = predict(S_smooth, noise_ratio) #pmin(turn_val, 1/(1-noise_ratio))
  final_S[noise_ratio >= turn_arg] = turn_val
  
  
  ## we are using sqrt(final_S) as the correction_factor here, but we refer to final_S as the correction factor in the paper.
  return(list(noise_ratio = noise_ratio, correction_factor = sqrt(final_S)))
}


