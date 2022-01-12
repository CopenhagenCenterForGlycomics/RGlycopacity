
test <- function() {
  allgenes = unique(glycogenome_ordering$HGNC)
  allgenes = sample(allgenes,30)
  #ggplot(data.frame(gene=allgenes,baz=sample(c('a','b'),length(allgenes),T),x=rnorm(length(allgenes),10,0.1)),aes(x,y=gene))+geom_point()+facet_glycogenome(genes=vars(gene))
  ggplot(data.frame(gene=rep(allgenes,3),baz=sample(c('a','b'),3*length(allgenes),T),x=rnorm(3*length(allgenes),10,0.1)),aes(x))+geom_histogram()+facet_glycogenome(genes=vars(gene))
}

rle_strings = function(strings) {
  rl <- rle(strings)
  lst <- lapply(split(seq_along(strings), rep(seq_along(rl$values), rl$lengths)), range)
  names(lst) <- rl$values
  do.call(rbind, lst)
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
         panel.spacing.y = unit(0.1,"lines"),
         strip.background.x = element_rect(size=1),
         strip.placement = "outside",
         strip.text.x = element_text(color = "black", face = "bold"),
         strip.background.y = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank()
    )

    params$rows = list(glycogene_pathway=quo(glycogene_pathway))

    drawn = ggplot2::FacetGrid$draw_panels(panels, layout, x_scales, y_scales, ranges, coord, data, theme, params)
    label_widths = lapply(drawn$grobs[which(grepl('axis-l',drawn$layout$name))], grid::grobWidth)

    used_pathways = unique(layout[,c('ROW','glycogene_pathway')])[,2]

    pathways_indices = rle_strings(stringr::str_replace(as.character(layout$glycogene_pathway),":.*",""))



    drawn$grobs[grepl('strip-l',drawn$layout$name)] = lapply(used_pathways, function(x) generate_pathway_element(x,extra_width=label_widths[[1]]))
    drawn$layout[grepl('strip-l',drawn$layout$name),]$clip = 'off'

    drawn = gtable::gtable_add_cols(drawn, unit(1, "cm"), 0)

    apply(cbind(pathways_indices,rownames(pathways_indices)), 1, function(run) {
      label = grid::textGrob(run[3],rot=90)
      line = grid::linesGrob(x = unit(c(0.75, 0.75), "npc"),
          y = unit(c(0.01, 0.99), "npc"),gp = grid::gpar(color = "#00000000",lwd=2))
      title = grid::grobTree(label,line)
      drawn <<- gtable::gtable_add_grob(drawn,grobs = title,t =as.numeric(run[1])*2,b=as.numeric(run[2])*2,l=1,r=1)
    })

    return(drawn)
  }

)

}