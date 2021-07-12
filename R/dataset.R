#' DPAGT1 normalised pseudobulks of glycogenes
#'
#' A dataset containing normalised pseudobulks of mouse and human PanglaoDB Chromium10x clusters for 224 glycogenes.
#'
#' @format A data frame with 969526 rows and 5 variables
#' \describe{
#'   \item{rel_diff}{DPAGT1 normalized pseudobulks (-7.627605--5.298152)}
#'   \item{tissue}{Tissue and cell type of each cluster, separated by colon (nTissue = 190; nCelltype = 101)}
#'   \item{gene}{The glycogene for which a pseudobulk is computed (nGene = 224)}
#'   \item{source}{Organism from which the cluster is derived (Human, Mouse)}
#'   \item{samples}{I have no idea at all}
#' }
#' @source repo:glyco_atlas_data/generate_log_dpagt1_norm; revision:07c39ba-dirty
"tissue_results_dpagt1"
#> [1] "tissue_results_dpagt1"
