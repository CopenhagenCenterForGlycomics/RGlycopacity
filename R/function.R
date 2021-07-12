#'@title shoud fit on one line, be written in sentence case, and not end in a full stop.
#'
#'@description this comes first in the documentation and should briefly describe what the function does
#'
#'@details this is a (often long) section that comes after the argument description and should provide any other important details of how the function operates. The details are optional.
#'
#'@return
#'
#'@examples
get_tissue_results_dpagt1 <- function() {
  # TODO load single cell results from a txt file or variable somewhere [tissue_results_dpagt1]
}

#'@title shoud fit on one line, be written in sentence case, and not end in a full stop.
#'
#'@description this comes first in the documentation and should briefly describe what the function does
#'
#'@details this is a (often long) section that comes after the argument description and should provide any other important details of how the function operates. The details are optional.
#'
#'@return
#'
#'@examples
get_glycogenes <- function() {
  # TODO pull glycogenes of interest from a txt file or variable somewhere
}

#'@title shoud fit on one line, be written in sentence case, and not end in a full stop.
#'
#'@description this comes first in the documentation and should briefly describe what the function does
#'
#'@details this is a (often long) section that comes after the argument description and should provide any other important details of how the function operates. The details are optional.
#'
#'@param data description
#'@param reference_data description
#'@param tissues description
#'@param quantile_cutoff description
#'@param reference_gene description
#'@param genes description
#'@param log_transform description
#'@param reference description
#'
#'@return
#'
#'@examples
#'
#'\dontrun{
#'encloses code that should not be run.
#'}
#'\dontshow{
#'encloses code that is invisible on help pages,
#'but will be run both by the package checking tools,
#'and the example() function.
#'This was previously testonly, and that form is still accepted.
#'}
#'\donttest{
#'encloses code that typically should be run,
#'but not during package checking.
#'The default run.donttest = interactive() leads example() use in other help page examples to skip  \donttest sections appropriately.
#'}
main <- function(data, reference_data = get_tissue_results_dpagt1, tissues = c(), quantile_cutoff = c(0.25, 0.75), reference_gene = c('DPAGT1'), genes = get_glycogenes, log_transform = c(log, log1p)) {

  cl <- match.call(expand.dots = T)
  cl_args <- names(cl)
  fx_args <- names(as.list(args(main)))
  match(cl_args, fx_args)
  for(arg in as.list(args(main))) { cl[[arg]] <- match.arg(arg) }

  ## retrieving genes to normalize ----
  genes <-
    switch(
      c('function','character')[c(
        is.function(eval(cl$genes)),
        is.character(cl$genes))],
      'function' = match.fun(cl$genes)(),
      'character' = unique(cl$genes, cl$reference_gene),
      rownames(data))

  ## ensuring genes of single cell reference and supplied samples are in order ----
  ordering <-
    sort(unique(c(
      rownames(data),
      rownames(tissue_results_dpagt1))))
  tissue_results_dpagt1 <- tissue_results_dpagt1[ordering,]
  data <- data[ordering,]

  ## filtering samples and reference by retrieved genes ----
  tissue_results_dpagt1 <- tissue_results_dpagt1[genes,]
  data <- data[genes,]

  ## filtering reference by choice of tissue ----
  # TODO pseudocode currently
  tissue_results_dpagt1 <-
    switch(
      c('subset')[length(tissues) > 0],
      'subset' = tissue_results_dpagt1[tissue_results_dpagt1$tissue %in% tissue,],
      tissue_results_dpagt1)

  ## computing value for specified quantile cutoffs on dpagt1 normalized data
  # TODO pseudocode currently
  expression_quantiles <- aggregate.data.frame(
    x = tissue_results_dpagt1,
    by = tissue_results_dpagt1$gene,
    FUN = ...,
    simplify = F)
  names(all_ranges) <- as.character(quantile_cutoff)

  ## computing log of reference gene in all samples ----
  reference <- switch(
    c('one','greater')[c(
      length(cl$reference_gene) == 1,
      length(cl$reference_gene) > 1)],
    'one' = data[reference_gene,],
    'greater' = exp(log(colMeans(data[reference_gene,]))))
  log_reference_data <- match.fun(cl$log_transform)(reference)

  ## computing log of all genes in all samples ---
  log_data <- match.fun(cl$log_transform)(data)

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
    normalised_expression,
    expression_quantiles,
    predicted_expression)

}
