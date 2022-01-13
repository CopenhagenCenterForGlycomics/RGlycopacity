generate_package_data = function() {
  glycogenome_ordering = read.delim('data/glycogenome_ordering.tsv',header=T)
  pathway_ordering = read.delim('data/pathway_ordering.tsv',header=T)
  glycogenome_ordering$Group.Name = factor(glycogenome_ordering$Group.Name, pathway_ordering$Group.Name)
  pathway_colours = read.delim('data/pathway_colours.tsv',header=T,comment.char='',na.strings=c(''))
  usethis::use_data(glycogenome_ordering,overwrite=T)
  usethis::use_data(pathway_colours,overwrite=T)
}