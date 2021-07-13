#'@title Retrieval of reference gene normalised pseudobulks from PanglaoDB
#'
#'@description this function retrieves the pre-computed reference gene normalised pseudobulks from PanglaoDB
#'
#'@param organism a character argument specifying the human or mouse specific reference dataset to return. Default  value is 'human'
#'
#'@return a data.frame of reference gene normalized pseudobulks
#'
#'@examples
#'\dontrun{
#'reference_data <- get_tissue_results_dpagt1()
#'reference_data <- get_tissue_results_dpagt1('mouse')
#'}
get_tissue_results_dpagt1 <- function(organism = c('human','mouse')) {

  organism <- match.arg(organism, choices = c('human','mouse'))

  reference_data <-
    switch(
      organism,
      'human' = subset(tissue_results_dpagt1, source == 'human'),
      'mouse' = subset(tissue_results_dpagt1, source == 'mouse'))

  reference_data

}

#'@title Retrieval of glycogenes from set list.
#'
#'@description this function retrieves the set of glycogenes to subset our reference and input data by
#'
#'@param organism a character argument specifying the human or mouse specific gene symbols to return. Default  value is 'human'
#'
#'@return a chqracter vector of glycogenes
#'
#'@examples
#'\dontrun{
#'genes_human <- get_glycogenes()
#'genes_mouse <- get_glycogenes('mouse')
#'}
get_glycogenes <- function(organism = c('human','mouse')) {

  organism <- match.arg(organism, choices = c('human','mouse'))

  switch(
    organism,
    'human' = human_glycogenes,
    'mouse' = mouse_glycogenes)

}

#'@title Predicting glycosylation capacity in bulk data from single cell pseudobulks
#'
#'@description this function determines whether a gene is expressed in bulk data by leveraging single cell reference data from PanglaoDB
#'
#'@param data an n gene x m sample matrix of bulk data values. Expression units can be in either TPM or FPKM
#'@param reference_data one of either a) a function returning a data.frame of reference gene normalized single cell pseudobulks, or b) a user supplied data.frame in the same format as that returned by the function in a. Default is (a)
#'@param organism an (a) character argument specifying the human or mouse specific datasets and gene symbols to return, or (b) an NA value indicating that user is using own data. Default value is 'human'
#'@param tissues one of either a) a character vector of tissues and/or celltypes that the reference data can be filtered for, or b) an NA value indicating that all tissues/celltypes should be used. Default is (b)
#'@param quantile_cutoff a numeric vector of values between (0, 1) indicating the quantile to be computed for the normalised pseudobulk distribution of each gene. Default is c(0.25, 0.75)
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
#'input_data <- as.matrix(read.table(
#'system.file("extdata", "GT_FPKM.csv", package = "glycoCapacityPredictR", mustWork = TRUE),
#'header = TRUE, sep = ',', row.names = 1))
#'
#'head(main(input_data)$normalised_expression)
#'
#'head(main(input_data,
#'quantile_cutoff = c(0.95))$normalised_expression)
#'
#'head(main(input_data,
#'reference_gene = c("DPAGT1","STT3A"))$normalised_expression)
#'
#'head(main(input_data,
#'tissues = 'Peripheral_blood_mononuclear_cells:Dendritic cells')$normalised_expression)
#'
#'head(main(input_data,
#'genes = c("XYLT1","XYLT2","STT3B"))$normalised_expression)
#'
#'head(main(input_data,
#'log_transform = 'log1p')$normalised_expression)
#'
#'@importFrom stats aggregate.data.frame quantile setNames
#'
#'@export
main <- function(data, reference_data = get_tissue_results_dpagt1, organism = c('human', 'mouse'), tissues = NA, quantile_cutoff = c(0.25, 0.75), reference_gene = c('DPAGT1'), genes = get_glycogenes, log_transform = c('log', 'log1p')) {

  ## stop condition for empty arguements ----
  null_args <- names(Filter(Negate(isFALSE), eapply(environment(), is.null)))
  if(length(null_args) > 0) {
    stop(paste('The following arguments require an input value:', paste0(null_args, collapse = ','))) }

  ## stop condition for data ----
  if(!is.matrix(data)) { stop('data must be an object of class matrix') }

  ## stop condition for reference_data ----
  if(sum(!is.data.frame(reference_data), !is.function(reference_data)) == 2) {
    stop('reference_data must be an object of class function or data.frame') }

  ## stop condition for organism ----
  organism <- match.arg(organism)
  if(sum(!is.character(organism), !is.na(organism)) == 2) {
    stop('organism must be an object of class character or NA') }

  ## stop condition for tissues ----
  if(sum(!is.character(tissues), !is.na(tissues)) == 2 ) {
    stop('tissues must be an object of class character or NA') }

  ## stop condition for quantile_cutoff ----
  if(!is.numeric(quantile_cutoff)) {
    stop('quantile_cutoff must be an object of class numeric') }

  ## stop condition for reference_gene ----
  if(!is.character(reference_gene)) {
    stop('reference_gene must be an object of class character') }

  ## stop condition for genes ----
  suppressWarnings(
    if(sum(!is.function(genes), !is.character(genes), !is.na(genes)) == 3) {
      stop('genes must be an object of class function, character, or NA') })

  ## stop condition for log_transform ----
  log_transform <- match.arg(log_transform)
  log_transform <- match.fun(log_transform)

  ## retrieving genes to normalize ----
  genes <-
    switch(
      c('function','character','NA')[c(
        is.function(genes),
        is.character(genes),
        suppressWarnings(is.na(genes)))],
      'function' = match.fun(genes)(organism),
      'character' = unique(c(genes, reference_gene)),
      'NA' = rownames(data))

  ## retrieving reference dataset ----
  reference_data <-
    switch(
      c('function','matrix')[c(
        is.function(reference_data),
        is.matrix(reference_data))],
      'function' = match.fun(reference_data)(organism),
      'matrix' = reference_data)

  ## filtering datasets for retrieved genes ----
  reference_data <- reference_data[reference_data$gene %in% genes,]
  data <- data[rownames(data) %in% genes,]

  ## filtering reference dataset for select tissues ----
  reference_data <-
    switch(
      c('subset','all')[c(
        is.character(tissues),
        is.na(tissues))],
      'subset' = reference_data[reference_data$tissue %in% tissues,],
      'all' = reference_data)

  ## computing value for specified quantile cutoffs on dpagt1 normalized data
  expression_quantiles <- aggregate.data.frame(
    x = reference_data$rel_diff,
    by = list(gene = reference_data$gene),
    FUN = function(x) { sapply(quantile_cutoff, quantile, x = x) },
    simplify = T)
  rownames(expression_quantiles) <- expression_quantiles$gene
  expression_quantiles$gene <- NULL
  expression_quantiles <- as.matrix(expression_quantiles)
  colnames(expression_quantiles) <- as.character(quantile_cutoff)

  ## ensuring genes of single cell reference and supplied samples are in order ----
  ordering <-
    sort(unique(c(
      rownames(data),
      rownames(expression_quantiles))))

  missing_data <- ordering[-which(ordering %in% rownames(data))]
  if(length(missing_data) > 0) {
    empty_data <- matrix(NA, nrow = length(ordering), ncol = ncol(data), dimnames = list(ordering, colnames(data)))
    empty_data[match(rownames(data), ordering),] <- data
    data <- empty_data
    remove(empty_data)
  }

  missing_quantiles <- ordering[-which(ordering %in% rownames(expression_quantiles))]
  if(length(missing_quantiles) > 0) {
    empty_quantiles <- matrix(NA, nrow = length(ordering), ncol = ncol(expression_quantiles), dimnames = list(ordering, colnames(expression_quantiles)))
    empty_quantiles[match(rownames(expression_quantiles), ordering),] <- expression_quantiles
    expression_quantiles <- empty_quantiles
    remove(empty_quantiles)
  }

  data <- data[ordering,]
  expression_quantiles <- expression_quantiles[ordering,]

  ## computing log of reference gene in all samples ----
  reference <- switch(
    c('one','greater')[c(
      length(reference_gene) == 1,
      length(reference_gene) > 1)],
    'one' = data[reference_gene,],
    'greater' = exp(log(colMeans(data[reference_gene,]))))
  log_reference_data <- log_transform(reference)

  ## computing log of all genes in all samples ---
  log_data <- log_transform(data)

  ## computing normalized expression ----
  normalised_expression <- log_data - log_reference_data

  ## computing predicted expression for each gene at each cutoff ----
  predicted_expression <-
    setNames(
      lapply(data.frame(normalised_expression), function(sample_library) {
        sample_prediction <-
          setNames(
            do.call(cbind.data.frame, lapply(
              data.frame(expression_quantiles), function(cutoffs) {
                ifelse(sample_library > cutoffs, 1, 0)
              })),
            colnames(expression_quantiles))
        rownames(sample_prediction) <- ordering
        sample_prediction
      }),
      colnames(normalised_expression))

  list(
    normalised_expression = normalised_expression,
    expression_quantiles = expression_quantiles,
    predicted_expression = predicted_expression)

}

# TODO choice of any reference gene will require normalization on log transformed data in data-raw/DATASET.R
# TODO log transformation of all genes will need to occur in data-raw/DATASET.R
# TODO enforce selection of only human or mouse subsets, no mixing or matching allowed
#
