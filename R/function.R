#'@title Retrieval of reference gene normalised pseudobulks from PanglaoDB
#'
#'@description this function retrieves the pre-computed reference gene normalised pseudobulks from PanglaoDB
#'
#'@return a data.frame of reference gene normalized pseudobulks
#'
#'@examples
#'reference_data <- get_tissue_results_dpagt1
get_tissue_results_dpagt1 <- function() { tissue_results_dpagt1 }

#'@title Retrieval of glycogenes from set list.
#'
#'@description this function retrieves the set of glycogenes to subset our reference and input data by
#'
#'@return a chqracter vector of glycogenes
#'
#'@examples
#'genes <- get_glycogenes()
get_glycogenes <- function() { glycogenes }

#'@title Predicting glycosylation capacity in bulk data from single cell pseudobulks
#'
#'@description this function determines whether a gene is expressed in bulk data by leveraging single cell reference data from PanglaoDB
#'
#'@param data an n gene x m sample matrix of bulk data values. Expression units can be in either TPM or FPKM
#'@param reference_data one of either a) a function returning a data.frame of reference gene normalized single cell pseudobulks, or b) a user supplied data.frame in the same format as that returned by the function in a. Default is (a)
#'@param tissues one of either a) a character vector of tissues and/or celltypes that the reference data can be filtered for, or b) an NA value indicating that all tissues/celltypes should be used. Default is (b)
#'@param quantile_cutoff a numeric vector of values between [0, 1] indicating the quantile to be computed for the normalised pseudobulk distribution of each gene. Default is c(0.25, 0.75)
#'@param reference_gene a character vector of reference gene(s) by which to normalise the input data matrix. If more than one reference gene is supplied, the geometric mean of normalized pseudobulks for the genes are computed in each sample. Default is c("DPAGT1")
#'@param genes one of either a) a function returning a character vector of glycogenes to be normalised, b) a user supplied character vector of genes, or c) an NA value indicating that all genes in the input data should be normalised. Default is (a)
#'@param log_transform a function of a) natural log, b) log1p. Additional log transform functions can be supplied by the user. Default is (a).
#'
#'@return a named list of three data.frames.
#'normalised_expression are the reference gene normalised values in each bulk input sample.
#'expression_quantiles are the set of values computed at each specified quantile on the normalised pseudobulk distribution of each gene in the reference data.
#'predicted_expression are the set of prediction statuses computed from applying the expression quantiles from the reference data to each bulk input sample.
#'
#'@examples
#'input_data <- ...
#'reference_data <- ...
#'main(input_data)
#'main(input_data, quantile_cutoff = c(0.95))
#'main(input_data, reference_gene = c("DPAGT1","STT3A"))
#'main(input_data, tissues = 'Brain')
#'main(input_data, reference_data = reference_data)
#'main(input_data, genes = c("XYLT1","XYLT2","STT3B"))
#'main(input_data, log_transform = log10)
#'
#'@export
main <- function(data, reference_data = get_tissue_results_dpagt1, tissues = NA, quantile_cutoff = c(0.25, 0.75), reference_gene = c('DPAGT1'), genes = get_glycogenes, log_transform = c(log, log1p)) {

  ## stop condition for empty arguements ----
  null_args <- names(Filter(Negate(isFALSE), eapply(environment(), is.null)))
  if(length(null_args) > 0) {
    stop(paste('The following arguments require an input value:', paste0(null_args, collapse = ','))) }

  ## stop condition for data ----
  if(!is.matrix(data)) { stop('data must be an object of class matrix') }

  ## stop condition for reference_data ----
  if(sum(!is.data.frame(reference_data), !is.function(reference_data) == 2)) {
    stop('reference_data must be an object of class function or data.frame') }

  ## stop condition for tissues ----
  if(!is.character(tissues)) {
    stop('tissues must be an object of class character') }

  ## stop condition for quantile_cutoff ----
  if(!is.numeric(quantile_cutoff)) {
    stop('quantile_cutoff must be an object of class numeric')}

  ## stop condition for reference_gene ----
  if(!is.character(reference_gene)) {
    stop('reference_gene must be an object of class character') }

  ## stop condition for genes ----
  if(sum(!is.function(genes), !is.character(genes), !is.na(genes)) == 3) {
    stop('genes must be an object of class function, character, or NA') }

  ## stop condition for log_transform ----
  log_transform <- match.arg(log_transform)

  ## retrieving genes to normalize ----
  genes <-
    switch(
      c('function','character','NA')[c(
        is.function(genes),
        is.character(genes),
        is.na(genes))],
      'function' = match.fun(genes)(),
      'character' = unique(genes, reference_gene),
      'NA' = rownames(data))

  ## retrieving reference dataset ----
  reference_data <-
    switch(
      c('function','matrix')[c(
        is.function(reference_data),
        is.matrix(reference_data))],
      'function' = match.fun(reference_data)(),
      'matrix' = reference_data)

  ## filtering datasets for retrieved genes ----
  reference_data <- reference_data[reference_data$genes %in% genes,]
  data <- data[genes,]

  ## filtering reference dataset for select tissues ----
  reference_data <-
    switch(
      c('subset')[length(tissues) > 0],
      'subset' = reference_data[reference_data$tissue %in% tissue,],
      reference_data)

  ## computing value for specified quantile cutoffs on dpagt1 normalized data
  expression_quantiles <- aggregate.data.frame(
    x = reference_data$rel_diff,
    by = list(gene = reference_data$gene),
    FUN = function(x) { sapply(quantile_cutoff, quantile, x = x) },
    simplify = T)
  rownames(expression_quantiles) <- expression_quantiles$gene
  expression_quantiles$gene <- NULL
  expression_quantiles <- as.matrix(expression_quantiles)
  colnames(all_ranges) <- as.character(quantile_cutoff)

  ## ensuring genes of single cell reference and supplied samples are in order ----
  ordering <-
    sort(unique(c(
      rownames(data),
      rownames(expression_quantiles))))
  expression_quantiles <- expression_quantiles[ordering,]
  data <- data[ordering,]

  ## computing log of reference gene in all samples ----
  reference <- switch(
    c('one','greater')[c(
      length(reference_gene) == 1,
      length(reference_gene) > 1)],
    'one' = data[reference_gene,],
    'greater' = exp(log(colMeans(data[reference_gene,]))))
  log_reference_data <- match.fun(log_transform)(reference)

  ## computing log of all genes in all samples ---
  log_data <- match.fun(log_transform)(data)

  ## computing normalized expression ----
  normalised_expression <- log_data - log_reference_data

  ## computing predicted expression for each gene at each cutoff ----
  predicted_expression <-
    setNames(
      lapply(normalised_expression, function(sample_library) {
        setNames(
          do.call(cbind.data.frame, lapply(
            expression_quantiles, function(cutoffs) {
              ifelse(sample_library > cutoffs, 1, 0)
              })),
          names(expression_quantiles))
        }),
      names(normalised_expression))

  list(
    normalised_expression = normalised_expression,
    expression_quantiles = expression_quantiles,
    predicted_expression = predicted_expression)

}
