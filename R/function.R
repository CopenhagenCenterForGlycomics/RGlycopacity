#'@title trimmed mean pseudobulk and pseudopresence calculation
#'
#'@description computing pseudobulks on seurat object
#'
#'@param obj.seurat seurat object with raw counts and cell types as identities
#'
#'@return a list of pseudobulks, pseudopresences and cluster sizes
#'
#'@examples
#'
#'@importFrom Matrix t sparseVector rowMeans rowSums
#'
#'@export
calculatePseudobulks <- function(obj.seurat,renormalize=FALSE) {

  if (!requireNamespace('Seurat',quietly=TRUE)) {
    stop('Seurat package is not installed')
  }

  if (!requireNamespace('pbmcapply',quietly=TRUE)) {
    stop('pbmcapply package is not installed')
  }

  requireNamespace('Matrix')

  

  # 1e4 and log1p transforming raw counts
  if (renormalize) {
    obj.seurat = Seurat::NormalizeData(obj.seurat, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')
  }
  obj.sm = Seurat::GetAssayData(obj.seurat, slot = 'data')
  obj.sm = Matrix::t(obj.sm)

  # retaining quants in 10th-90th percentile per gene
  col_indices <- 1:ncol(obj.sm)
  row_indices <- 1:nrow(obj.sm)

  batches.list <- split(col_indices, ceiling(seq_along(col_indices)/1000))

  message("Batching Pseudobulk operation into",length(batches.list),"batches")

  pbmcapply_function <- function(i, c.df, r.vec) {

    value.vec <- c.df[r.vec, i]
    quant.vec <- value.vec[which(value.vec > 0)]
    cutof.vec <- quantile(quant.vec, na.rm = T, probs = c(0.1,0.9))
    cutof.max <- which(value.vec > 0 & (value.vec < cutof.vec[1] | value.vec > cutof.vec[2]))
    value.vec[cutof.max] <- NA
    value.vec <- Matrix::sparseVector(value.vec, seq(value.vec), length(value.vec))
    return(value.vec)

  }

  new.sm <-
    do.call(
      cbind,
      lapply(batches.list, function(batch) {
        do.call(
          cbind,
          lapply(
            pbmcapply::pbmclapply(
              X = batch,
              FUN = pbmcapply_function,
              obj.sm,
              row_indices,
              ignore.interactive = T),
            as, "sparseMatrix"))
        }))

  rownames(new.sm) <- rownames(obj.sm)
  colnames(new.sm) <- colnames(obj.sm)
  obj.sm <- Matrix::t(new.sm)

  # retrieving unique clusters
  cluster_ids <- as.character(unique(SeuratObject::Idents(obj.seurat)))

  # computing pseudobulk per cluster
  pseudobulks.sm <-
    SeuratObject::as.sparse(
      do.call(
        cbind,
        lapply(cluster_ids, function(cluster) {
          Matrix::rowMeans(obj.sm[, SeuratObject::WhichCells(obj.seurat, ident = cluster), drop = F], na.rm = T) })))
  colnames(pseudobulks.sm) <- cluster_ids

  # computing pseudopresence per cluster
  pseudopresence.sm <-
    SeuratObject::as.sparse(
      do.call(
        cbind,
        lapply(cluster_ids, function(cluster) {
          Matrix::rowSums(obj.sm[, SeuratObject::WhichCells(obj.seurat, ident = cluster), drop = F] > 0, na.rm = T) })))
  colnames(pseudopresence.sm) <- cluster_ids

  # retrieving cluster sizes
  cluster_sizes.vec <- table(SeuratObject::Idents(obj.seurat))[cluster_ids]

  list(
    pseudobulks.sm = pseudobulks.sm,
    pseudopresence.sm = pseudopresence.sm,
    cluster_sizes.vec =  cluster_sizes.vec)

}

#'@title single cell gene expression prediction
#'
#'@description predicted expression of a gene in a single cell cluster
#'
#'@param ncells number of cells in the cluster
#'@param pseudobulk mean expression of the gene in the cluster
#'
#'@return a pvalue corresponding to the probability that the gene is expressed in the cluster
#'
#'@examples
#'
#'@importFrom stats predict pt
expression_prediction <- function(ncells, pseudobulk) {

  s <-
    stats::predict(
      expression_model.lm$model2.lm,
      newdata = data.frame(x = ncells),
      se.fit = TRUE,
      level = 0.95,
      interval = "confidence")$fit[,'fit']

  n   <- 7 # chosen based on the number of pseudobulk sizes that model is fit to
  ts  <- (mean(expression_model.lm$xy.df$y[1:7]) - pseudobulk)/(s/sqrt(n))
  p   <- 2*pt(-abs(ts), n-1)

  return(p)

}

#'@export
predictExpressionOverCutoff <- function(single_cell.list) {
  single_cell.list$pseudobulks.sm > mean(expression_model.lm$xy.df$y[1:7])
}

#'@title Make a data frame predicting the expression of glycogenes from pseudobulk data
#'
#'@description Using the result data from calculatePseudobulks, we can
#'             calculate the set of genes that are predicted to be expressed, and
#'             where their pseudobulk values sit compared to a reference range
#'@export
calculateExpressedRange <- function(single_cell.list,reference=available_ranges()) {
#  What would the exp do here?
#  normalised = normaliseOnSEGs(exp(single_cell.list$pseudobulks.sm+1),log_transform=F)
  normalised = normaliseOnSEGs(single_cell.list$pseudobulks.sm,log_transform=F)
  intersects = calculateNormalisedIntersects(normalised,reference)
  intersects = ifelse(as.matrix(single_cell.list$pseudobulks.sm[rownames(intersects),]) > mean(RGlycopacity:::expression_model.lm$xy.df$y[1:7]), intersects,NA)
  as.data.frame(intersects)
}

#'@title single cell gene expression prediction wrapper
#'
#'@description predicted expression of a gene in a single cell cluster
#'
#'@param single_cell.list list of pseudobulk sparse matrix and vector of cluster sizes
#'
#'@return an n gene x k cluster matrix of pvalues corresponding to the probability that the gene is expressed in the cluster
#'
#'@examples
#'
#'@export
predictMinimalExpression <- function(single_cell.list) {

  expression.mat <-
    do.call(
      cbind,
      mapply(
        expression_prediction,
        as.list(single_cell.list$cluster_sizes.vec),
        data.frame(single_cell.list$pseudobulks.sm, check.rows = F),
        SIMPLIFY = F))

  rownames(expression.mat) <- rownames(single_cell.list$pseudobulks.sm)
  colnames(expression.mat) <- colnames(single_cell.list$pseudobulks.sm)

  expression.mat

}

#'@title single cell pseudopresence prediction
#'
#'@description predicted expression of pseudobulk stability in a single cell cluster
#'
#'@param pseudobulk mean expression of the gene in the cluster
#'@param pseudopresence number of cells with at least one transcript for the gene in the cluster
#'
#'@return a list of expected mean pseudopresence, expected var pseudopresence and a pvalue corresponding to the probability that the pseudobulk is stable
#'
#'@examples
#'
#'@importFrom stats predict pt
pseudopresence_prediction <- function(pseudobulk, pseudopresence = NULL) {

  mu_model.lm <- pseudopresence_model.lm$glycogene_mean_fit.lm
  sd_model.lm <- pseudopresence_model.lm$glycogene_sd_fit.lm

  # computing sample population parameter mean ----
  mu <-
    stats::predict(
      mu_model.lm,
      newdata = data.frame(pseudobulk_mean = pseudobulk),
      se.fit = TRUE,
      level = 0.95,
      interval = "confidence")$fit[,'fit']

  # computing sample population parameter variance ----
  s <-
    stats::predict(
      sd_model.lm,
      newdata = data.frame(pseudobulk = pseudobulk),
      se.fit = TRUE,
      level = 0.95,
      interval = "confidence")$fit[,'fit']

  if(!is.null(pseudopresence)) {

    # n chosen based on the number of pseudobulks that model is fit to
    n <- nrow(pseudopresence_model.lm$glycogene_sd_fit.lm$model)
    ts  <- (mu - pseudopresence)/(s/sqrt(n))
    p   <- 2*pt(-abs(ts), n-1)

  } else { p <- NA }

  list(
    expected.mu = mu,
    expected.s = s,
    pval = p)
}

#'@title single cell pseudopresence prediction wrapper
#'
#'@description predicted expression of pseudobulk stability in a single cell cluster
#'
#'@param single_cell.list list of pseudobulk and pseudopresence sparse matrices
#'
#'@return a list of matrices corresponding to  expected mean pseudopresence, expected var pseudopresence and stability pvalue
#'
#'@examples
#'
#'@export
predictPseudopresenceStability <- function(single_cell.list) {

  pseudopresence.list <-
    mapply(
      pseudopresence_prediction,
      data.frame(single_cell.list$pseudobulks.sm),
      data.frame(single_cell.list$pseudopresence.sm),
      SIMPLIFY = F)

  expected_mu.mat <- do.call(cbind, lapply(pseudopresence.list, `[[`, 1))
  rownames(expected_mu.mat) <- rownames(single_cell.list$pseudobulks.sm)
  colnames(expected_mu.mat) <- names(pseudopresence.list)

  expected_var.mat <- do.call(cbind, lapply(pseudopresence.list, `[[`, 2))
  rownames(expected_var.mat) <- rownames(single_cell.list$pseudobulks.sm)
  colnames(expected_var.mat) <- names(pseudopresence.list)

  pval.mat <- do.call(cbind, lapply(pseudopresence.list, `[[`, 3))
  rownames(pval.mat) <- rownames(single_cell.list$pseudobulks.sm)
  colnames(pval.mat) <- names(pseudopresence.list)

  pseudopresence.list <- list(
    expected_mu.mat = expected_mu.mat,
    expected_var.mat = expected_var.mat,
    pval.mat = pval.mat)

  pseudopresence.list

}

#'@title reference gene normalisation
#'
#'@description normalising on a reference gene for downstream comparison of glycogene ranges of expression
#'
#'@param genes a character vector of the sample genes
#'@param values a named numeric vector of the counts to be transformed. Names are genes
#'@param log_transform boolean indicating if a natural log + 1 transform should be applied
#'
#'@return a data.frame of genes, annotations, and normalised counts
#'
seg_normalisation <- function(
  genes,
  values,
  log_transform = T) {

  if(log_transform) { values <- log(values + 1) }

  seg_indices <- which(genes %in% segs)

  seg_values <- values[seg_indices]

  geom_mean <- if(any(seg_values < 0)) { mean() } else { function(x) { exp(mean(log(x))) } }

  seg_value <- geom_mean(seg_values)

  values <- values - seg_value

  values
}

#'@title reference gene normalisation wrapper
#'
#'@description normalising on a reference gene for downstream comparison of glycogene ranges of expression
#'
#'@param input.mat n gene x k sample matrix of suitable counts
#'@param genes a named numeric vector of the counts to be transformed. default is rownames of input_mat
#'@param log_transform boolean indicating if a natural log + 1 transform should be applied
#'
#'@return a data.frame of genes, annotations, and normalised counts
#'
#'@examples
#'
#'@export
normaliseOnSEGs <- function(input.mat, genes = row.names(input.mat), log_transform = T) {

  rel_diff.mat <-
    do.call(
      cbind,
      lapply(
        data.frame(as.matrix(input.mat)),
        seg_normalisation,
        genes = genes,
        log_transform = log_transform))

  rownames(rel_diff.mat) <- rownames(input.mat)
  rownames(rel_diff.mat) <- rownames(input.mat)

  rel_diff.mat
}

#'@title compute clrs for a vector of values
#'
#'@description this function determines the clrs of a vector
#'
#'@param values a vector of values
#'
#'@return a vector of clrs
#'
#'@examples
#'
#'@importFrom compositions clr
clr_transform <- function(values) { compositions::clr(values) }

#'@title wrapper for computing clrs on test data
#'
#'@description this function determines the clrs of a gene in test data
#'
#'@param input.mat an n gene x k sample matrix of values
#'
#'@return the set of ranges computed from applying the specified quantiles to each gene in the test data
#'
#'@examples
#'
#'@export
transformToCLR <- function(input.mat) {

  clrs.mat <-
    as.matrix(
      do.call(
        rbind,
        lapply(
          data.frame(t(input.mat)),
          clr_transform)))

  colnames(clrs.mat) <- colnames(input.mat)
  rownames(clrs.mat) <- rownames(input.mat)

  clrs.mat
}

#'@title compute quantiles for a vector of values
#'
#'@description this function determines the quantiles of a vector
#'
#'@param values a vector of values
#'@param quantiles a vector of desired quantiles
#'
#'@return a vector of quantiles
#'
#'@examples
#'
range_generation <- function(values, quantiles) { quantile(x = values, probs = quantiles, type = 1, na.rm = T) }

#'@title Predicting dynamic ranges in test data
#'
#'@description this function determines the ranges of a gene in test data
#'
#'@param input.mat an n gene x m sample matrix of values
#'@param quantiles desired ranges to compute
#'
#'@return the set of ranges computed from applying the specified quantiles to each gene in the test data
#'
#'@examples
#'
#'@export
getNormalisedRanges <- function(input.mat, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)) {

  ranges.mat <-
    as.matrix(
      do.call(
        rbind,
        lapply(
          data.frame(t(input.mat)),
          range_generation,
          quantiles = quantiles)))

  colnames(ranges.mat) <- quantiles
  rownames(ranges.mat) <- rownames(input.mat)

  ranges.mat
}

#'@title Available reference ranges
#'@export
available_ranges = function() {
c('bulk_ranges',
        'tabula_sapiens_ranges',
        'panglao_musculus_ranges',
        'TCGA_ranges',
        'GTEx_ranges',
        'panglao_musculus_unknown_ranges',
        'panglao_musculus_known_ranges')
}

#'@title Predicting glycosylation capacity in test data from reference data
#'
#'@description this function determines whether a gene is expressed in test data by leveraging reference data
#'
#'@param input.mat an n gene x m sample matrix of values.
#'@param dynamic_ranges.char character reference to a stored set of ranges
#'
#'@return the set of prediction statuses computed from applying the expression quantiles from the reference data to each sample in the test data
#'
#'@examples
#'
#'@importFrom stats quantile setNames
#'@importFrom Matrix t rowSums
#'
#'@export
calculateNormalisedIntersects <- function(input.mat, dynamic_ranges.char = available_ranges()) {

  ## stop condition for empty arguements ----
  null_args <- names(Filter(Negate(isFALSE), eapply(environment(), is.null)))
  if(length(null_args) > 0) {
    stop(paste('The following arguments require an input value:', paste0(null_args, collapse = ','))) }

  ## stop condition for test_data ----
  if(!is.matrix(input.mat)) { stop('test_data must be an object of class matrix') }

  dynamic_ranges.char <-
    match.arg(
      dynamic_ranges.char
    )

  ## filtering datasets for retrieved genes ----
  shared_genes <- intersect(unique(rownames(input.mat)), unique(colnames(dynamic_ranges[[dynamic_ranges.char]]$expression_quantiles)))

  ## ensuring genes of reference_data and test_data are in order ----
  test_data <- input.mat[rownames(input.mat) %in% shared_genes, , drop = F]
  reference_data <- Matrix::t(dynamic_ranges[[dynamic_ranges.char]]$expression_quantiles[, shared_genes, drop = F])
  reference_data <- as.matrix(reference_data)
  reference_data <- reference_data[rownames(test_data), , drop = F]

  ## computing predicted expression for each gene at each cutoff ----
  predicted_expression <-
    mapply(
      function(sample_library, cutoffs) { Matrix::rowSums(sample_library > cutoffs) },
      data.frame(test_data),
      replicate(ncol(test_data), reference_data, simplify = F))

  ## labeling each prediction ----
  predicted_classification <-
    matrix(
      as.character(
        factor(
          predicted_expression,
          levels = 0:11,
          labels = c('<0%', paste0('>',seq(0,100,10))))),
      nrow = nrow(test_data), ncol = ncol(test_data),
      dimnames = dimnames(test_data), byrow = F)

  predicted_classification

}
