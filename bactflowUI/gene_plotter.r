#!/bin/env Rscript

crans = c('optparse', 'tidyverse', 'readr', 'ape', 'BiocManager', 'pheatmap', 'glue')

for(pkg in crans){
    if(!(requireNamespace(pkg, quietly = TRUE))){
        options(repos = c(CRAN = "https://cloud.r-project.org/"))
        install.packages(pkg, quiet = TRUE, keep_outputs = F)
    }
}


if(!(requireNamespace('MatrixGenerics', quietly = T))){
    BiocManager::install('MatrixGenerics', force = T)
}

for (pkgs in c(crans, 'MatrixGenerics')){
    library(pkgs, character.only = TRUE, verbose = FALSE)
}

# library(optparse)
# library(tidyverse)
# library(readr)
# library(ape)
# library(MatrixGenerics)
# library(pheatmap)
# library(glue)

# gene classification

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

g_count = read_tsv(file = opt$count_table) %>% column_to_rownames("Product")
products = rownames(g_count)

if(!dir.exists(opt$output_dir)){
    dir.create(opt$output_dir)
}

enzymes = read_tsv(opt$enzymes)



all = enzymes$Enzymes




# mat1 = apply(g_count, 2, function(x){log(1+x)})
# mat1 <- g_count
# mat1 <- as.matrix(mat1)

# sds <- MatrixGenerics::rowSds(mat1)
# o <- order(sds, decreasing = TRUE)
# dist_mat_gene = dist(mat1[o,], method = "euclidean")
# dist_mat_genome = dist(t(mat1[o,]), method = "euclidean")

# h_1 <- hclust(dist_mat_gene, method = "ward.D2")#for column in this case
# h_2 <- hclust(dist_mat_genome, method = "ward.D2")





# heat_p =  pheatmap(mat = t(mat1[o,]), angle_col = 45, fontsize_row = 12,
#                       border_color = NA, 
#                       cellwidth = 22, cellheight = 20, cutree_cols = 4, cutree_rows = 4, fontsize_number = 8,
#                      cluster_row = h_2, 
#                      cluster_cols = h_1,
#                       col = RColorBrewer::brewer.pal(9, "Reds"), 
#                       main = "Heat map of log-scaled (all)genes abundance of different genomes")


# ggsave(plot = heat_p, filename = paste0(opt$output_dir, "all_genes_abundance.jpeg"), device = "jpeg", dpi = 300, height =20, width = 50, limitsize = F)




### Now for the genes of interst



disc_color = function(n, method = "brewer", seed = 1990){
    if(method == "rainbow"){
        return(rainbow(n))
    } else if(method == "brewer"){
        
        
        if(n<=11){
            cols = RColorBrewer::brewer.pal(n = 11, name = "Set3")[1:n]
          

             
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
            cols = cols[permute::shuffle(cols)][1:n]

        } else if(!is.numeric(n)){
            stop(
                "Please use a valid number!"
            )
        }
    } else if (method == "viridis"){
        cols = viridis::viridis(n = n, option = "D")
    }
return(cols)
}

matches <- list()
orig_name <- list()
for (i in seq_along(all)) {
  found <- unique(products[grep(all[i], products, ignore.case = TRUE)])  
  if (length(found) > 0) {
    matches[[i]] <- found  
    orig_name[[i]] <- all[i]
  } else {
    matches[[i]] <- NA  
    orig_name[[i]] <- NA
  }
}


# making col annot

matches = unlist(matches)
matches = matches[!is.na(matches)] %>% unique()
orig_name = unlist(orig_name)
orig_name = orig_name[!is.na(orig_name)]

wanted_enz = enzymes


wanted_enz = wanted_enz[wanted_enz$Enzymes %in% orig_name,]
col_annot = data.frame(enz = matches, Functions = NA, stringsAsFactors = F)



for(i in seq_along(matches)){
 mf = unique(wanted_enz$Function[sapply(wanted_enz$Enzymes, function(x) grepl(x, col_annot$enz[i], ignore.case = TRUE))])
 col_annot$Functions[i] <-  ifelse(length(mf) > 0, paste(mf, collapse = "|"), NA)
}


col_annot = col_annot %>% distinct(enz, Functions)  %>% column_to_rownames('enz')


# mat1 = apply(plantarum, 2, function(x){log(1+x)})
mat1 <- g_count
mat1 <- as.matrix(mat1)


mat_gen = mat1[matches,] 
outname_table = "table_of_gene_abundance.tsv"
outname_jpeg = "requested_genes_abundance.png"
title_p = glue("Heat map of  {length(matches)} genes abundance in different {opt$name_organism} genomes")

if(opt$prevalence == TRUE){
    mat_gen[mat_gen>0] <- 1 
    

}

#  mat1 = scale(mat1) # optional
sds <- MatrixGenerics::rowSds(mat_gen)
o <- order(sds, decreasing = TRUE)
dist_mat_gene = dist(mat_gen[o,], method = "euclidean")
dist_mat_genome = dist(t(mat_gen[o,]), method = "euclidean")

h_1 <- hclust(dist_mat_gene, method = "ward.D2")#for column in this case
h_2 <- hclust(dist_mat_genome, method = "ward.D2")


cols = matrix( disc_color( n=2*length(unique(col_annot[["Functions"]])))[sort(seq(length(unique(col_annot[["Functions"]]))), decreasing = T)], dimnames = list(unique(col_annot$Functions)))

if(opt$prevalence == TRUE){
    outname_table = "table_of_gene_prevalence.tsv"
    outname_jpeg = "requested_genes_prevalence.png"
    title_p = glue("Heat map of {length(matches)} genes prevalence (contingency) in {opt$name_organism} genomes")

    heat_p =  pheatmap(mat = t(mat_gen[o,]), angle_col = 45, fontsize_row = 12, 
                      border_color = "#fffefd", col = c("#b5b5b5", "#009275"),
                      annotation_col = col_annot,
                      annotation_colors = list(
                        Functions = cols[,1]
                      ),
                      legend_breaks = seq(from = range(mat_gen)[1], to = range(mat_gen)[2], 1),
                      cellwidth = 25, cellheight = 20, cutree_cols = 4, cutree_rows = 3, fontsize_number = 8,
                    # cluster_cols = F, 
                    # cluster_rows = F,
                     cluster_row = h_2, 
                     display_numbers = t(mat_gen[o,]),
                     cluster_cols = h_1,
                     number_color  =  '#ffffff',
                      # col = RColorBrewer::brewer.pal(9, "Blues"), 
                      # col = disc_color(n = range(mat_gen)[2]+1),
                      main = title_p, silent = TRUE)

} else {
    heat_p =  pheatmap(mat = t(mat_gen[o,]), angle_col = 45, fontsize_row = 12, 
                      annotation_col = col_annot,
                      annotation_colors = list(
                        Functions = cols[,1]
                      ),
                      legend_breaks = seq(from = range(mat_gen)[1], to = range(mat_gen)[2], 1),
                      cellwidth = 25, cellheight = 20, cutree_cols = 4, cutree_rows = 3, fontsize_number = 8,
                     cluster_row = h_2, 
                     display_numbers = t(mat_gen[o,]),
                     cluster_cols = h_1,
                     number_color  =  '#090909',
                      # col = RColorBrewer::brewer.pal(9, "Blues"), 
                      col = disc_color(n = range(mat_gen)[2]+1),
                      main = title_p, silent = TRUE)
}





write.table(x = mat_gen[o,] , file  = glue("{opt$output_dir}/{outname_table}"), sep  = '\t', row.names = T)



ggsave(plot = heat_p, filename = glue("{opt$output_dir}/{outname_jpeg}"), device = "png", dpi = 300, height =opt$height, width = opt$width, limitsize = F)
