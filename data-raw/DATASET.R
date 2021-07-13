## code to prepare `tissue_results_dpagt1` dataset goes here

# Log-DPAGT1 data gen

## Log-DPAGT1-Normalised data generation

# LOAD ----
library(dplyr)
library(compositions)

# FUNCTIONS ----

## expression model and cutoff used for computing threshold and pvalue ----
expression_model.lm <- readRDS('lm_chromium10x_mupb_norm_nziqr_threshold_cutoff.Rds')
cutoff.int <- mean(expression_model.lm$xy.df$y[1:7])

## computing threshold on each cluster ----
compute_threshold <- function(ncells) {

  s <-
    stats::predict(
      expression_model.lm$model2.lm,
      newdata = data.frame(x = ncells),
      se.fit = TRUE,
      level = 0.95,
      interval = "confidence")$fit[,'fit']

  return(s)

}

## computing pvalue for each glycogene on each cluster
compute_pvalue <- function(pseudobulk, threshold) {

  n   <- 7 # chosen based on the number of pseudobulk sizes that model is fit to
  ts  <- (cutoff.int - pseudobulk)/(threshold/sqrt(n))
  p   <- 2*pt(-abs(ts), n-1)

  return(p)

}

## wrappers for compute_pvalue and compute_threshold ----
expression_status_call <- function(threshold,pseudobulk) {
  if (is.na(pseudobulk)) {
    return("unexpressed")
  }
  if (pseudobulk == 0) {
    return("unexpressed")
  }

  expression_status_call_raw(threshold,pseudobulk)
}

## wrappers for compute_pvalue and compute_threshold ----
expression_status_call_raw <- function(threshold, pseudobulk) {

  expression <- pseudobulk > cutoff.int
  pvalue     <- compute_pvalue(pseudobulk, threshold)
  if (is.na(expression)) {
    return('unexpressed')
  }
  status <- if(expression) { 'expressed' } else { 'unexpressed' }
  status <- if(pvalue > 5E-03) { 'unpredictable' } else { status }

  return(status)
}

## wrappers for compute_pvalue and compute_threshold ----
is_expressed <- function(pseudobulks) {
  cells = names(pseudobulks)
  expressed = sapply(cells, function(cell) { expression_status_call( attributes(alldata)$thresholds[cell], pseudobulks[cell] ) })
  ifelse(expressed == 'expressed', TRUE, ifelse(expressed == 'unexpressed',FALSE,NA))
}


# MAIN ----

## reading in mouse and human gene symbols ----
ortho_trans_table=with(read.delim('mouse_orthology_human_genes.tsv'), setNames(mouse_symbol,human_symbol))
glycogenes=read.delim('glycogenes.tsv',header=F)$V1

## loading mouse datasets ----
alldata = readRDS('alldata_mouse.Rds')
alldata_clean_mouse = alldata[!grepl('SMALL',rownames(alldata)),]
attributes(alldata_clean_mouse)$clustersize = attributes(alldata)$clustersize
colnames(alldata_clean_mouse) = recode(colnames(alldata_clean_mouse), Tmem5='Rxylt1',Tmem246='Pgap4',Kdelc1='Poglut2',Kdelc2='Poglut3')
rm(alldata)

## storing mouse clustersize and threshold as attributes in empty list ----
alldata=list()
attributes(alldata)$clustersize = c(attributes(alldata_clean_mouse)$clustersize)
attributes(alldata)$thresholds = sapply(attributes(alldata)$clustersize,compute_threshold)

## filtering for clusters with Dpagt1 and Alg2 and Alg13 or Alg14 expression in mouse ----
required_genes=c('Dpagt1','Alg2')
passing_cells = apply(apply(alldata_clean_mouse[,required_genes],2,is_expressed) ,1,all)
passing_cells = passing_cells & apply(apply(alldata_clean_mouse[,c('Alg13','Alg14')],2,is_expressed),1,any)
alldata_clean_mouse = alldata_clean_mouse[which(passing_cells),]

## natural log transform of mouse pseudobulks ----
alldata_clean_mouse_log = apply(alldata_clean_mouse,2,log)
alldata_clean_mouse_log[!is.finite(alldata_clean_mouse_log)] <- 0

## TODO description
mouse_melt = subset(reshape2::melt(alldata_clean_mouse_log[,intersect(ortho_trans_table[glycogenes],colnames(alldata_clean_mouse_log))]),value != 0)
mouse_melt$Var2 = setNames(names(ortho_trans_table),ortho_trans_table)[as.character(mouse_melt$Var2)]
mouse_melt = cbind(subset(mouse_melt, Var2 != 'cell'),data='mouse')
rm(alldata_clean_mouse)
rm(alldata_clean_mouse_log)

## loading human datasets ----
alldata = readRDS('alldata_human.Rds')
alldata_clean_human = alldata[!grepl('SMALL',rownames(alldata)),]
attributes(alldata_clean_human)$clustersize = attributes(alldata)$clustersize
colnames(alldata_clean_human) = recode(colnames(alldata_clean_human), TMEM5='RXYLT1',TMEM246='PGAP4',KDELC1='POGLUT2',KDELC2='POGLUT3')
rm(alldata)

## storing human clustersize and threshold as attributes in empty list ----
alldata=list()
attributes(alldata)$clustersize = c(attributes(alldata_clean_human)$clustersize)
attributes(alldata)$thresholds = sapply(attributes(alldata)$clustersize,compute_threshold)

## filtering for clusters with DPAGT1 and ALG2 and ALG13 or ALG14 expression in human ----
required_genes=c('DPAGT1','ALG2')
passing_cells = apply(apply(alldata_clean_human[,required_genes],2,is_expressed) ,1,all)
passing_cells = passing_cells & apply(apply(alldata_clean_human[,c('ALG13','ALG14')],2,is_expressed),1,any)
alldata_clean_human = alldata_clean_human[which(passing_cells),]

## natural log transform of human pseudobulks ----
alldata_clean_human_log = apply(alldata_clean_human,2,log)
alldata_clean_human_log[!is.finite(alldata_clean_human_log)] <- 0

## TODO description
human_melt = cbind(subset(reshape2::melt(alldata_clean_human_log[,intersect(glycogenes,colnames(alldata_clean_human_log))]),value != 0),data='human')
human_melt$Var2 = as.character(human_melt$Var2)
rm(alldata_clean_human)
rm(alldata_clean_human_log)

## combining natural log transformed human and mouse data ----
full_log = rbind(human_melt,mouse_melt)
full_log$tissue_cell = stringr::str_replace(stringr::str_replace(full_log$Var1,'^[^:]+:',''),'_cluster.*','')
full_log$cell = stringr::str_replace( full_log$tissue_cell,'.*:','')
rm(human_melt)
rm(mouse_melt)

## normalising natural log transformed pseudobulks to DPAGT1 ----
log_dpagt1_norm = merge(data.frame(Var1=full_log[full_log$Var2 == 'DPAGT1',]$Var1, dpagt=full_log[full_log$Var2 == 'DPAGT1',]$value),full_log,by='Var1')
log_dpagt1_norm$value = log_dpagt1_norm$value - log_dpagt1_norm$dpagt
log_dpagt1_norm = setNames(log_dpagt1_norm[,c('value','tissue_cell','Var2','data')],c('rel_diff','tissue','gene','source'))
log_dpagt1_norm$samples = 100

## normalising natural log transformed cutoff to DPAGT1
attributes(log_dpagt1_norm)$cutoffs = with(subset(full_log,Var2=='DPAGT1'),setNames(log(cutoff.int)-value,tissue_cell))

## additional attributes regarding repo
# attributes(log_dpagt1_norm)$repo = "glyco_atlas_data/generate_log_dpagt1_norm"
# attributes(log_dpagt1_norm)$revision = system("git describe --always --dirty",intern=T)
# saveRDS(log_dpagt1_norm,'pseudobulk_log_dpagt1_norm.Rds')

## renaming variable for shits and giggles ----
human_glycogenes <- sort(unique(subset(log_dpagt1_norm, source == 'human')$gene))
mouse_glycogenes <- sort(unique(subset(log_dpagt1_norm, source == 'mouse')$gene))
tissue_results_dpagt1 <- log_dpagt1_norm

## writing variables as internal package .rda
usethis::use_data(human_glycogenes, mouse_glycogenes, tissue_results_dpagt1, internal = T, overwrite = TRUE)
