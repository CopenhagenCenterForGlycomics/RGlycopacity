
test <- function() {
  ggplot(data.frame(gene=c('GALNT2','GALNT10','POFUT1','POFUT2'),baz=rep(c('a','b'),4),x=1:4),aes(x,y=gene))+geom_point()+facet_glycogenome(genes=vars(gene))
}

if (requireNamespace("ggplot2",quietly=TRUE)) {

facet_glycogenome <- function(genes = NULL, cols = NULL, scales = "free",
                       space = "free", shrink = TRUE
                       ) {

  scales <- match.arg(scales, c("fixed", "free_x", "free_y", "free"))
  free <- list(
    x = any(scales %in% c("free_x", "free")),
    y = any(scales %in% c("free_y", "free"))
  )

  space <- match.arg(space, c("fixed", "free_x", "free_y", "free"))
  space_free <- list(
    x = any(space %in% c("free_x", "free")),
    y = any(space %in% c("free_y", "free"))
  )

  facets_list = list(
    rows = ggplot2:::compact_facets(ggplot2:::as_facets_list(genes)),
    cols = ggplot2:::compact_facets(ggplot2:::as_facets_list(cols))
  )

  ggproto(NULL, FacetGlycogenome,
    shrink = shrink,
    params = list(rows = facets_list$rows, cols = facets_list$cols, margins = FALSE,
      free = free, space_free = space_free, labeller = "label_value",
      as.table = TRUE, switch = "both", drop = TRUE)
  )
}

FacetGlycogenome <- ggplot2::ggproto("FacetGlycogenome", ggplot2::FacetGrid,
  compute_layout = function(data,params) {
    data = lapply(data, function(df) {
      genes = combine_vars(data, params$plot_env, params$rows, drop = params$drop)
      df$glycogene_pathway = with(RGlycopacity::glycogenome_ordering, setNames(Group.Name,HGNC))[genes[[1]]]
      df
    })
    params$rows = list(glycogene_pathway=quo(glycogene_pathway))
    return(ggplot2::FacetGrid$compute_layout(data,params))
  },
  map_data = function(data,layout,params) {
    genes = combine_vars(list(data), params$plot_env, params$rows, drop = params$drop)
    data$glycogene_pathway = with(RGlycopacity::glycogenome_ordering, setNames(Group.Name,HGNC))[genes[[1]]]
    data$gene = factor(data$gene,rev( unique(RGlycopacity::glycogenome_ordering$HGNC )))
    params$rows = list(glycogene_pathway=quo(glycogene_pathway))
    return(ggplot2::FacetGrid$map_data(data,layout,params))
  },
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord, data, theme, params) {
    theme = theme + theme(panel.spacing.x = unit(0.1,"lines"),
         panel.spacing.y = unit(1,"lines"),
         strip.background.x = element_rect(size=1),
         strip.placement = "outside",
         strip.text.x = element_text(color = "black", face = "bold"),
         strip.background.y = element_blank(),
         strip.text.y = element_text(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank()
    )
    params$rows = list(glycogene_pathway=quo(glycogene_pathway))
    return(ggplot2::FacetGrid$draw_panels(panels, layout, x_scales, y_scales, ranges, coord, data, theme, params))
  }

)

}