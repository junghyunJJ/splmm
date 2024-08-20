# rm(list = ls())

# library(tidyverse)
# library(data.table)
# library(jjutil)

# # c("1-1", "2-3", "2-5", "2-8", "18-64", "T4857")
# # # control	1-1
# # # Mid-AD	2-3
# # # control	2-5
# # # Mid-AD	2-8
# # # control	18-64
# # # Mid-AD	T4857

# # sel_sample <- "1-1"
# # indir <- "data"

# # # phenotype (expression)
# # raw_exp <- fread(str_glue("{indir}/exp_{sel_sample}.cvs.gz")) %>% jjutil::convert_dataframe()
# # raw_exp %>% h
# # raw_exp %>% dim

# # exp <- raw_exp[apply(raw_exp, 1, sum) != 0, ] %>% t %>% scale %>% t
# # exp %>% dim
# # exp %>% h
# # # rowMeans(exp) %>% mean
# # # apply(exp, 1, mean)
# # # apply(exp, 1, sd)

# # # genetic effect -> spatial effect (diatance kernel)
# # distance_kernel <- fread(str_glue("{indir}/K_{sel_sample}.cvs.gz")) %>% jjutil::convert_dataframe()
# # all.equal(colnames(exp), rownames(distance_kernel))

# # distance_kernel %>% h
# # distance_kernel %>% dim

# # # Env effect -> cell type effect (cell type proportion from cell2location)
# # cell2loc <- fread(str_glue("data/{sel_sample}_Cell2location_results.csv.gz")) %>% jjutil::convert_dataframe()
# # cell2loc %>% head
# # rownames(cell2loc) <- sub(str_glue("{sel_sample}_"), "", rownames(cell2loc))
# # all.equal(rownames(cell2loc), rownames(distance_kernel))


# # cell2loc %>% head

# # celltype_kernel <- cell2loc %>% cal_linear_kernel
# # celltype_kernel %>% h
# # celltype_kernel %>% dim
# # celltype_kernel %>% head()

# indir <- "/home/jungj2/project/test/test_Celina/Real data analysis"
# # load in data
# load(str_glue("{indir}/scRNA_data_RCC_PD47171in_tumor_interface.RData"))
# # single cell reference data
# load(str_glue("{indir}/data_RCC_PD47171in_tumor_interface.RData"))
# # 10x visium spatial transcriptomics data
# names(sc_section)
# names(data_section)

# # filter genes: 
# # 1. remove mt genes:
# mito_genes <- unique(c(grep("^MT-", rownames(data_section$count)), grep("^mt-", rownames(data_section$count))))
# # 2. remove genes expressed in less than 20 locations
# low_genes <- which(rowSums(as.matrix(data_section$count) > 0) < 20)
# remove_genes <- unique(c(mito_genes, low_genes))
# data_expr_in <- data_section$count[-remove_genes, ]


# ################################################################
# #### 1. prep ###################################################
# ################################################################

# data_section$celltype_proportion %>% head
# data_section$location_coord %>% head


# Obj <- CELINA::Create_Celina_Object(
#   celltype_mat = t(data_section$celltype_proportion),
#   gene_expression_mat = as.matrix(data_expr_in),
#   location = as.matrix(data_section$location_coord),
#   covariates = NULL,
#   project = "Kidney Cancer"
# )

# # if you want to further filter cell types based on their total proportion across spots, or you only want to test a subset of cell types, you can select the cell types here:
# filtered_cell_types <- colnames(data_section$celltype_proportion)[which(colSums(data_section$celltype_proportion) > (dim(data_section$celltype_proportion)[1] * 0.01))]
# filtered_cell_types

# Obj <- CELINA::preprocess_input(Obj,
#   # Celina object
#   cell_types_to_test = filtered_cell_types,
#   # a vector of cell types to be used for testing
#   scRNA_count = as.matrix(sc_section$count),
#   # a gene x cell expression matrix of reference scRNA-seq data
#   sc_cell_type_labels = sc_section$meta$broad_type,
#   # a vector of cell type labels for each cell in scRNA_count
#   threshold = 5e-5
# )


# celltype_mat <- Obj@celltype_mat %>% t
# celltype_mat %>% dim
# celltype_mat %>% head

# exp <- as.matrix(Obj@gene_expression_mat)
# exp %>% h
# exp %>% dim
# all.equal(rownames(celltype_mat), colnames(exp))

# # location <- as.matrix(Obj@location)
# location <- as.matrix(data_section$location_coord) #NOTE!!! original coord data
# all.equal(rownames(celltype_mat), rownames(location))
# save.image("tmp.RData")
# ################################################################
# #### 2. Celina #################################################
# ################################################################

# # Obj <- CELINA::Calculate_Kernel(Obj, approximation = FALSE)
# # Obj %>% glimpse
# # Obj <- CELINA::Testing_interaction_all(Obj, celltype_to_test = "RCC", num_cores = 40)


################################################################
#### 2. run ####################################################
################################################################
# load("tmp.RData")
# saveRDS(exp, "data/example_exp.rds")
# saveRDS(location, "data/example_location.rds")
# saveRDS(celltype_mat, "data/example_celltype_mat.rds")

rm(list = ls())

library(tidyverse)
library(data.table)
library(jjutil)

exp <- readRDS("data/example_exp.rds")
location <- readRDS("data/example_location.rds")
celltype_mat <- readRDS("data/example_celltype_mat.rds")
exp %>% dim


source("R/splmm.R")
source("R/kernel.R")
source("R/utils.R")

sel_sample <- "kidneycancer"
sel_celltype <- "RCC"
exp %>% dim


res <- splmm(exp = exp, coord = location, celltype_prop = celltype_mat, 
  sel_celltype = "RCC", 
  sel_gene = rownames(exp)[1],
  # sel_gene = c("FABP7", "GC", "SLC13A1", "CXCL9", "PAX8", "RPL28"),
  bandwidthtype = "Scott",
  path_mtg = "./mtg2",
  tmpdir = str_glue("tmp_{sel_sample}"),
  nthread = 20, verbose = 1, remove_tmpdir = FALSE
)

getwd()

################################################################
#### 3. old run ################################################
################################################################

celltype_kernel <- celltype_mat %>% cal_linear_kernel()

sel_celltype_kernel <- celltype_mat[, colnames(celltype_mat) %in% sel_celltype, drop = FALSE] %>% cal_linear_kernel()
nonsel_celltype_kernel <- celltype_mat[, !colnames(celltype_mat) %in% sel_celltype, drop = FALSE] %>% cal_linear_kernel()

# j <- 11
# splmm(unlist(exp[j, ]), distance_kernel, celltype_kernel,
#   path_mtg = "./mtg2",
#   tmpdir = str_glue("tmp/tmp_{sel_sample}_{j}"),
#   nthread = 1, verbose = 1, remove_tmpdir = TRUE
# )

# j <- 1
# splmm(unlist(exp[j, ]), distance_kernel, sel_celltype_kernel,
#   path_mtg = "./mtg2",
#   tmpdir = str_glue("tmp/tmp_{sel_sample}_{j}"),
#   nthread = 1, verbose = 1, remove_tmpdir = TRUE
# )

# j <- 1
# splmm(unlist(exp[j, ]), distance_kernel, nonsel_celltype_kernel,
#   path_mtg = "./mtg2",
#   tmpdir = str_glue("tmp/tmp_{sel_sample}_{j}"),
#   nthread = 1, verbose = 1, remove_tmpdir = TRUE
# )
