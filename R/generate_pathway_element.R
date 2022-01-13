
generate_pathway_element <- function(label,extra_width=grid::unit(0,"npc")) {


  colours = RGlycopacity::pathway_colours

  bar_lookup = with(colours,setNames(bar,pathway))
  shade_lookup = with(colours,setNames(shade,pathway))
  shape_lookup = with(colours,setNames(sugar,pathway))

  pathway_class = stringr::str_replace(label,":.*","")
  pathway = stringr::str_replace(label,".*:","")
  if (is.na(pathway)) {
    return(grid::nullGrob())
  }

  shade = NULL
  if (pathway_class != "Initiation") {
    shade = shade_lookup[pathway]
  }
  pathway_strip = grid::nullGrob()
  if (!is.na(bar_lookup[pathway])) {
    pathway_strip = grid::rectGrob(grid::unit(0,"npc"), grid::unit(0.5,"npc"), grid::unit(0.5,"npc"), grid::unit(1,"npc"), just=c(0,0.5), gp = grid::gpar(fill=bar_lookup[pathway],lwd=1), name="path_strip")
  }

  min_shade_width = unit(0.5,"cm")

  pathway_shade = grid::rectGrob(grid::unit(0.5,"npc"), grid::unit(0.5,"npc"), min_shade_width + extra_width , grid::unit(1,"npc"), just=c(0,0.5), gp = grid::gpar(fill=shade,lwd=0), name="path_shade")

  sugar_grob = grid::nullGrob()
  if (requireNamespace('ggsugar',quietly=TRUE)) {
    if (!is.na(shape_lookup[pathway])) {
      sug = ggsugar:::get_template_sugar(shape_lookup[pathway])
      sugar_viewport = grid::viewport(
        x=grid::unit(0.5,"npc") + grid::unit(0,"native"),
        y=grid::unit(0,"native"),
        width=min_shade_width,
        height=grid::unit(1,"npc"),
        just=c(0,0)
      )
      sugar_grob = grid::gTree(vp=sugar_viewport,children = do.call(grid::gList,sug))
    }
  }

  grid::grobTree(pathway_shade,pathway_strip,sugar_grob)
}