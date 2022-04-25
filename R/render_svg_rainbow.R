v8_ctx = NULL

get_v8 = function() {
  if (is.null(v8_ctx)) {
    v8_ctx <- V8::v8(console=T)
    v8_ctx$source(system.file("sviewer-headless.js", package = "ggsugar", mustWork = TRUE));
    v8_ctx$source(system.file("glycotopiary.js", package = "RGlycopacity", mustWork = TRUE));
  }
  v8_ctx
}

svg_doc_portrait=paste(readLines(system.file("glycotopiary_portrait.svg", package = "RGlycopacity", mustWork = TRUE)),collapse='\n')
svg_doc_landscape=paste(readLines(system.file("glycotopiary_landscape.svg", package = "RGlycopacity", mustWork = TRUE)),collapse='\n')
reactions_json=paste(readLines(system.file("reactions.json", package = "RGlycopacity", mustWork = TRUE)),collapse='\n')

perform_glycotopiary = function(svg_doc, reactions_json,wanted_genes=NULL) {
  v8_ctx = get_v8();
  v8_ctx$assign("svg_doc",svg_doc);
  v8_ctx$assign("reactions",reactions_json);
  v8_ctx$assign("wanted_genes",wanted_genes);
  v8_ctx$eval(paste("perform_glycotopiary(svg_doc,reactions,wanted_genes).then( res => console.r.assign('svg_temp',res) )",sep=""));
  retval=get('svg_temp',envir=.GlobalEnv);
  rm('svg_temp',envir=.GlobalEnv);
  retval
}

recolour_labels = function(svg_doc, colours) {
  v8_ctx = get_v8();
  v8_ctx$assign("svg_doc",svg_doc);
  v8_ctx$assign("colours",colours);
  v8_ctx$eval(paste("recolour_labels(svg_doc,colours).then( res => console.r.assign('svg_temp',res) )",sep=""));
  retval=get('svg_temp',envir=.GlobalEnv);
  rm('svg_temp',envir=.GlobalEnv);
  retval
}

set_size = function(svg_doc, portrait=TRUE) {
  v8_ctx = get_v8();
  v8_ctx$assign("svg_doc",svg_doc);
  v8_ctx$assign("portrait",portrait);
  v8_ctx$eval(paste("set_size(svg_doc,portrait).then( res => console.r.assign('svg_temp',res) )",sep=""));
  retval=get('svg_temp',envir=.GlobalEnv);
  rm('svg_temp',envir=.GlobalEnv);
  retval
}

render_svg = function(svg_doc) {
  input_svg = tempfile("input",fileext=".svg")
  cat(svg_doc,file=input_svg,sep="\n")
  cairo_svg = tempfile("cairo",fileext=".svg")
  suppressMessages(grConvert::convertPicture(input_svg,cairo_svg))
  grImport2::readPicture(cairo_svg)
  grImport2::pictureGrob(grImport2::readPicture(cairo_svg))
}