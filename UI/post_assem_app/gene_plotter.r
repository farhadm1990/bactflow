#!/bin/env Rscript

needed = c(
  'optparse', 'readr', 'dplyr', 'tidyr', 'tibble', 'ggplot2',
  'stringr', 'purrr', 'ape', 'pheatmap', 'glue', 'RColorBrewer'
)

for (pkg in needed) {
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
        install.packages(pkg, quiet = TRUE)
    }
}

if (!(requireNamespace('MatrixGenerics', quietly = TRUE))) {
    if (!(requireNamespace('BiocManager', quietly = TRUE))) {
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
        install.packages('BiocManager', quiet = TRUE)
    }
    BiocManager::install('MatrixGenerics', force = TRUE, ask = FALSE, update = FALSE)
}

for (pkg in c(needed, 'MatrixGenerics')) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}


options_list <- list(
    make_option(c("-c", "--count_table"), type = 'character', default = "./count_mixed.tsv", help = "Path to the count matrix of annotated genes [default: %default]"),
    make_option(c("-e", "--enzymes"), type = "character", default = "wanted_enzymes.tsv", help = "A tab seperated table with two column, Function and Enzymes. [default: %default]"),
    make_option(c('-p', "--prevalence"), type="logical", default=FALSE, help = "Convert gene abundance into prevalence with only presence 1/absence 0.. [default: %default]"),
    make_option(c('-o', '--output_dir'), type = 'character', default = getwd()),
    make_option(c('-n', '--name_organism'), type = "character", default = "Bacteria", help =  "Species name, e.g. Plantarum"),
    make_option(c('-w', '--width'), type ="double", default=20, help = "Width of the output plot in cm [default: %default]"),
    make_option(c('-l', '--height'), type ="double", default=20, help = "Height of the output plot in cm [default: %default]")
)

opt <- parse_args(OptionParser(option_list = options_list))

norm_str <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- tolower(x)
  x <- gsub("\u03b2|ß", "beta", x, perl = TRUE)
  x <- gsub("\u03b1", "alpha", x, perl = TRUE)
  x <- gsub("\u03b3", "gamma", x, perl = TRUE)
  x <- gsub("[\u2212\u2013\u2014]", "-", x, perl = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

# Enzyme query must appear inside Bakta Product/Gene (fixed, normalized).
enzyme_in <- function(enzyme, field) {
  q <- norm_str(enzyme)
  t <- norm_str(field)
  if (nchar(q) < 3 || !nzchar(t)) return(FALSE)
  grepl(q, t, fixed = TRUE)
}

# Gene symbol may appear inside the enzyme label (e.g. ruvB in "... subunit RuvB").
gene_in_enzyme <- function(gene, enzyme) {
  g <- norm_str(gene)
  q <- norm_str(enzyme)
  if (nchar(g) < 3 || !nzchar(q)) return(FALSE)
  grepl(g, q, fixed = TRUE)
}

g_count = read_tsv(file = opt$count_table, show_col_types = FALSE) %>% column_to_rownames("Product")
products = rownames(g_count)

if(!dir.exists(opt$output_dir)){
    dir.create(opt$output_dir, recursive = TRUE)
}

enzymes = read_tsv(opt$enzymes, show_col_types = FALSE)
if (!all(c("Function", "Enzymes") %in% names(enzymes))) {
  stop("Enzyme file must have columns 'Function' and 'Enzymes'.")
}
enzymes = enzymes %>%
  mutate(Enzymes = as.character(Enzymes), Function = as.character(Function)) %>%
  filter(!is.na(Enzymes), nzchar(Enzymes))

all = enzymes$Enzymes

# Optional Gene↔Product map written by gene_counter_bakta.py
map_path <- sub("\\.tsv$", ".gene_map.tsv", opt$count_table)
gene_map <- NULL
if (file.exists(map_path)) {
  gene_map <- read_tsv(map_path, show_col_types = FALSE) %>%
    mutate(
      Gene = ifelse(is.na(Gene), "", as.character(Gene)),
      Product = ifelse(is.na(Product), "", as.character(Product))
    )
  message("Loaded gene map: ", map_path, " (", nrow(gene_map), " rows)")
}

disc_color = function(n, method = "brewer", seed = 1990){
    if(method == "rainbow"){
        return(rainbow(n))
    } else if(method == "brewer"){
        if(n<=11){
            cols = RColorBrewer::brewer.pal(n = max(3, n), name = "Set3")[1:n]
        } else if( n>11 && n <=74){
            sets_name  =  RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% filter(category == "qual") %>% dplyr::select(names) %>% pull
            sets_num = RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% filter(category == "qual") %>% dplyr::select(maxcolors) %>% pull
            pooled = list()
            for(i in 1:length(sets_name)){
                pooled[[i]] <- RColorBrewer::brewer.pal(n = as.numeric(sets_num[i]), name = sets_name[i])
            }
            cols = do.call(c, pooled)[1:n]
        } else if( n > 74){
            sets_name  =  RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names") %>% dplyr::select(names) %>% pull
            sets_num = RColorBrewer::brewer.pal.info %>% data.frame() %>% rownames_to_column("names")  %>% dplyr::select(maxcolors) %>% pull
            pooled = list()
            for(i in 1:length(sets_name)){
                pooled[[i]] <- RColorBrewer::brewer.pal(n = as.numeric(sets_num[i]), name = sets_name[i])
            }
            cols = do.call(c, pooled)
            set.seed(seed)
            cols = cols[sample.int(length(cols), n, replace = TRUE)]
        } else if(!is.numeric(n)){
            stop("Please use a valid number!")
        }
    } else if (method == "viridis"){
        cols = viridis::viridis(n = n, option = "D")
    }
return(cols)
}

# Match each wanted enzyme against Bakta Product names and Gene symbols.
matches <- character(0)
orig_name <- character(0)
match_notes <- character(0)

for (i in seq_along(all)) {
  enz <- all[i]
  found_prod <- products[vapply(products, function(p) enzyme_in(enz, p), logical(1))]
  found_via_gene <- character(0)
  if (!is.null(gene_map) && nrow(gene_map) > 0) {
    keep <- vapply(seq_len(nrow(gene_map)), function(j) {
      g <- gene_map$Gene[j]
      p <- gene_map$Product[j]
      isTRUE(enzyme_in(enz, p)) ||
        isTRUE(enzyme_in(enz, g)) ||
        isTRUE(gene_in_enzyme(g, enz))
    }, logical(1))
    # Guard: never keep NA logicals (would select unintended rows)
    keep[is.na(keep)] <- FALSE
    hit_rows <- gene_map[keep, , drop = FALSE]
    found_via_gene <- unique(hit_rows$Product[nzchar(hit_rows$Product)])
  }
  found <- unique(c(found_prod, found_via_gene))
  found <- found[!is.na(found) & found %in% products]
  if (length(found) > 0) {
    matches <- c(matches, found)
    orig_name <- c(orig_name, rep(enz, length(found)))
    match_notes <- c(match_notes, paste0(enz, " -> ", paste(found, collapse = "; ")))
  }
}

matches = unique(matches)
orig_name = unique(orig_name[!is.na(orig_name)])

if (length(matches) == 0) {
  stop(
    "No enzymes from the wanted list were found in Bakta Product/Gene fields.\n",
    "Check that enzyme names resemble Bakta Product names or Gene symbols ",
    "(e.g. 'beta-glucosidase' or 'ruvB').\n",
    "Count table: ", opt$count_table, "\n",
    "Enzyme file: ", opt$enzymes
  )
}

message("Matched ", length(matches), " Bakta product(s) for ", length(orig_name), " enzyme quer(ies):")
for (note in unique(match_notes)) message("  ", note)

wanted_enz = enzymes[enzymes$Enzymes %in% orig_name, , drop = FALSE]
col_annot = data.frame(enz = matches, Functions = NA_character_, stringsAsFactors = FALSE)

for (i in seq_along(matches)) {
  mf = unique(wanted_enz$Function[vapply(wanted_enz$Enzymes, function(x) {
    isTRUE(enzyme_in(x, col_annot$enz[i]))
  }, logical(1))])
  col_annot$Functions[i] <- ifelse(length(mf) > 0, paste(mf, collapse = "|"), NA_character_)
}

col_annot = col_annot %>% distinct(enz, Functions) %>% column_to_rownames('enz')

mat1 <- as.matrix(g_count)
# Keep a true 2D matrix even when only one gene or one genome matches
mat_gen = mat1[matches, , drop = FALSE]
storage.mode(mat_gen) <- "numeric"

outname_table = "table_of_gene_abundance.tsv"
outname_jpeg = "requested_genes_abundance.png"
title_p = glue("Heat map of  {length(matches)} genes abundance in different {opt$name_organism} genomes")

if (isTRUE(opt$prevalence)) {
    mat_gen[mat_gen > 0] <- 1
    outname_table = "table_of_gene_prevalence.tsv"
    outname_jpeg = "requested_genes_prevalence.png"
    title_p = glue("Heat map of {length(matches)} genes prevalence (contingency) in {opt$name_organism} genomes")
}

n_genes <- nrow(mat_gen)
n_genomes <- ncol(mat_gen)

# Variance ranking (safe for 1-row / 1-col matrices)
if (n_genes >= 1 && n_genomes >= 1) {
  sds <- MatrixGenerics::rowSds(mat_gen, useNames = TRUE)
  sds[is.na(sds)] <- 0
  o <- order(sds, decreasing = TRUE)
} else {
  o <- seq_len(n_genes)
}

mat_ord <- mat_gen[o, , drop = FALSE]
annot_ord <- col_annot[rownames(mat_ord), , drop = FALSE]

cluster_genes <- n_genes >= 2
cluster_genomes <- n_genomes >= 2
h_1 <- if (cluster_genes) hclust(dist(mat_ord, method = "euclidean"), method = "ward.D2") else FALSE
h_2 <- if (cluster_genomes) hclust(dist(t(mat_ord), method = "euclidean"), method = "ward.D2") else FALSE

cut_cols <- if (cluster_genes) min(4, n_genes) else NA
cut_rows <- if (cluster_genomes) min(3, n_genomes) else NA

fn_levels <- unique(annot_ord$Functions[!is.na(annot_ord$Functions)])
if (length(fn_levels) == 0) {
  annot_ord$Functions <- "Unassigned"
  fn_levels <- "Unassigned"
}
fn_cols <- disc_color(n = max(3, length(fn_levels)))[seq_along(fn_levels)]
names(fn_cols) <- fn_levels
annot_colors <- list(Functions = fn_cols)

rng <- range(mat_ord, na.rm = TRUE)
legend_breaks <- if (isTRUE(opt$prevalence) || diff(rng) < 1) {
  sort(unique(as.numeric(mat_ord)))
} else {
  seq(from = rng[1], to = rng[2], by = 1)
}

heat_args <- list(
  mat = t(mat_ord),
  angle_col = 45,
  fontsize_row = 12,
  annotation_col = annot_ord,
  annotation_colors = annot_colors,
  cellwidth = 25,
  cellheight = 20,
  fontsize_number = 8,
  cluster_rows = h_2,
  cluster_cols = h_1,
  display_numbers = t(mat_ord),
  main = title_p,
  silent = TRUE
)

if (!is.na(cut_cols) && cut_cols >= 2) heat_args$cutree_cols <- cut_cols
if (!is.na(cut_rows) && cut_rows >= 2) heat_args$cutree_rows <- cut_rows

if (isTRUE(opt$prevalence)) {
  heat_args$border_color <- "#fffefd"
  heat_args$col <- c("#b5b5b5", "#009275")
  heat_args$number_color <- "#ffffff"
  heat_args$legend_breaks <- legend_breaks
} else {
  n_cols <- max(3, as.integer(rng[2]) + 1)
  heat_args$col <- disc_color(n = n_cols)
  heat_args$number_color <- "#090909"
  heat_args$legend_breaks <- legend_breaks
}

heat_p = do.call(pheatmap, heat_args)

write.table(x = mat_ord, file = glue("{opt$output_dir}/{outname_table}"), sep = '\t', row.names = TRUE)

ggsave(
  plot = heat_p,
  filename = glue("{opt$output_dir}/{outname_jpeg}"),
  device = "png",
  dpi = 300,
  height = opt$height,
  width = opt$width,
  limitsize = FALSE
)

message("Wrote ", glue("{opt$output_dir}/{outname_table}"), " and ", glue("{opt$output_dir}/{outname_jpeg}"))
