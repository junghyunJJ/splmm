rm(list = ls())

library(tidyverse)
library(data.table)
library(jjutil)

source("R/splmm.R")
# control	1-1
# Mid-AD	2-3
# control	2-5
# Mid-AD	2-8
# control	18-64
# Mid-AD	T4857

c("1-1", "2-3", "2-5", "2-8", "18-64", "T4857")

sel_sample <- "1-1"
indir <- "data"

# phenotype (expression)
raw_exp <- fread(str_glue("{indir}/exp_{sel_sample}.cvs.gz")) %>% jjutil::convert_dataframe()
raw_exp %>% h
raw_exp %>% dim

exp <- raw_exp[apply(raw_exp, 1, sum) != 0, ] %>% t %>% scale %>% t
exp %>% dim
exp %>% h
# rowMeans(exp) %>% mean
# apply(exp, 1, mean)
# apply(exp, 1, sd)

# genetic effect -> spatial effect (diatance kernel)
distance_kernel <- fread(str_glue("{indir}/K_{sel_sample}.cvs.gz")) %>% jjutil::convert_dataframe()
all.equal(colnames(exp), rownames(distance_kernel))

distance_kernel %>% h
distance_kernel %>% dim

# Env effect -> cell type effect (cell type proportion from cell2location)
cell2loc <- fread(str_glue("data/{sel_sample}_Cell2location_results.csv.gz")) %>% jjutil::convert_dataframe()
cell2loc %>% head
rownames(cell2loc) <- sub(str_glue("{sel_sample}_"), "", rownames(cell2loc))
all.equal(rownames(cell2loc), rownames(distance_kernel))

celltype_kernel <- cell2loc %>% cal_linear_kernel
celltype_kernel %>% h
celltype_kernel %>% dim

################################################################
#### 1. prep ###################################################
################################################################

# path_mtg <- "mtg2/mtg2.22_src/mtg2"
# tmpdir <- "./example"
# nthread <- 20
# verbose <- 1
# remove_tmpdir <- TRUE
# i <- 1
# res <- spLMM(unlist(exp[i, ]), distance_kernel, celltype_kernel, path_mtg = path_mtg, tmpdir = tmpdir, nthread = 20, verbose = 1, remove_tmpdir = TRUE)

# sub_exp <- exp[1:5, ]

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
idxs <- parallel::splitIndices(nrow(exp), 1000)[[i]]

# j <- 26
# j <- 11
res <- lapply(idxs, function(j) {
  splmm(unlist(exp[j, ]), distance_kernel, celltype_kernel,
        path_mtg = "mtg2/mtg2.22_src/mtg2",
        tmpdir = str_glue("tmp/tmp_{sel_sample}_{j}"),
        nthread = 1, verbose = 1, remove_tmpdir = TRUE
  )
})

if (!dir.exists(str_glue("res/splmm/{sel_sample}"))) {
  dir.create(str_glue("res/splmm/{sel_sample}"), recursive = TRUE)
}

saveRDS(res, str_glue("res/splmm/{sel_sample}/res_splmm_{sel_sample}_{i}.rds"))
