
#' @export
glycotopiaryGrob = function(expressed=glycogenome_ordering$HGNC,colours=NULL,portrait=TRUE) {
  if (requireNamespace('V8',quietly = TRUE) && requireNamespace('ggsugar',quietly = TRUE) && requireNamespace('grConvert',quietly = TRUE)) {
    result = perform_glycotopiary(ifelse(portrait,svg_doc_portrait, svg_doc_landscape),reactions_json,expressed)
    result = recolour_labels(result, colours)
    result = set_size(result, portrait)
    return( render_svg(result) )
  } else {
    stop("Method glycotopiaryGrob requires packages V8, ggsugar, grConvert and grImport2 to be installed")
  }
}

#' @export
glycotopiarySVG = function(expressed=glycogenome_ordering$HGNC,colours=NULL,portrait=TRUE) {
  if (requireNamespace('V8',quietly = TRUE) && requireNamespace('ggsugar',quietly = TRUE) ) {
    result = perform_glycotopiary(ifelse(portrait,svg_doc_portrait, svg_doc_landscape),reactions_json,expressed)
    result = recolour_labels(result, colours)
    result = set_size(result, portrait)
    return( result )
  } else {
    stop("Method glycotopiarySVG requires packages V8 and ggsugar to be installed")
    
  }
}