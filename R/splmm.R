library(tidyverse)
library(data.table)

std <- function(dat) {

  ncov <- ncol(dat)

  # cov <- scale(dat) %>% data.frame
  # standardization
  cov <- apply(dat, 2, function(x) {
    ave <- mean(x, na.rm = TRUE)
    std <- sd(x, na.rm = TRUE)
    (x - ave) / (sd(x) * sqrt(ncov))
    # (x - ave) / sd(x)
  })
  colnames(cov) <- NULL
  rownames(cov) <- NULL

  return(as.matrix(cov))
}


cal_linear_kernel <- function(dat) {
  stdcov <- std(dat)
  K <- tcrossprod(as.matrix(stdcov)) #K <- (as.matrix(stdcov) %*% t(as.matrix(stdcov)))
  # K <- (as.matrix(stdcov) %*% t(as.matrix(stdcov))) / ncol(stdcov)
  return(K)
}


# NOTE!!!! We need to update the function
cal_spatial_kernel <- function(normalized_expr, coord, kerneltype = "gaussian", bandwidthtype = "Silverman", bandwidth.set.by.user = NULL, sparseKernel = FALSE, sparseKernel_tol = 1e-20, sparseKernel_ncore = 1) {
  # cal spatial Kernel using SpatialPCA r package
  # https://lulushang.org/SpatialPCA_Tutorial/slideseq.html

  # normalized_expr: g x n (we used SCTransform)
  # coord: n x 2 (i.e., x and y)

  # The type of bandwidth to be used in Gaussian kernel,
  #   1. "SJ" for Sheather & Jones (1991) method (usually used in small size datasets),
  #   2. "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually used in large size datasets).

  # scale expr data to calculate "bandwidth"
  expr <- as.matrix(as.data.frame(t(scale(t(normalized_expr)))))

  if (is.null(bandwidth.set.by.user)) {
    bandwidth <- SpatialPCA::bandwidth_select(expr, method = bandwidthtype)
    cat("bandwidth: ", bandwidth, "\n", sep = "")

  } else {
    bandwidth <- bandwidth.set.by.user
    cat("bandwidth by user: ", bandwidth, "\n", sep = "")
  }

  # scale coordinate data
  location_normalized <- scale(coord)

  if (sparseKernel == FALSE) {
    kernelmat <- SpatialPCA::kernel_build(kerneltype = kerneltype, location = location_normalized, bandwidth = bandwidth)
  } else if (sparseKernel == TRUE) {
    kernelmat <- SpatialPCA::kernel_build_sparse(
      kerneltype = kerneltype,
      location = location_normalized, bandwidth = bandwidth,
      tol = sparseKernel_tol, ncores = sparseKernel_ncore
    )
  }
  return(kernelmat)
}


make_longform <- function(K) {
  n_sample <- nrow(K)
  rownames(K) <- seq(1, n_sample)
  colnames(K) <- seq(1, n_sample)

  K[upper.tri(K)] <- NA
  K <- as.data.frame(K)
  K <- cbind(id = rownames(K), K)
  K <- gather(K, "iid", "value", -id)
  K <- K %>%
    filter(!is.na(value)) %>%
    mutate(id = as.integer(id)) %>%
    mutate(iid = as.integer(iid)) %>%
    arrange(id, iid)

  return(K)
}


cal_Q <- function(K, A) {
  K <- as.matrix(K)
  A <- as.matrix(A)

  if (!matrixcalc::is.positive.definite(K)) {
    cat(" / Generate near PD: distance kernel", sep = "")
    K <- as.matrix(Matrix::nearPD(K)$mat)
  }
  K_chol <- t(chol(K))

  if (!matrixcalc::is.positive.definite(A)) {
    cat(" / Generate near PD: celltype kernel", sep = "")
    A <- as.matrix(Matrix::nearPD(A)$mat)
  }
  A_chol <- t(chol(A))

  Q <- (A_chol %*% t(K_chol)) + t(A_chol %*% t(K_chol))
  cat("\n")
  # if (!matrixcalc::is.positive.definite(Q)) {
  #     cat("Q: generate near PD...")
  #     Q <- as.matrix(Matrix::nearPD(Q)$mat)
  # }
  return(Q)
}


summry_res <- function(res) {
  # Vg
  res1 <- strsplit(res[7], ":") %>%
    unlist() %>%
    trimws()
  h2_1 <- sub("SE", "", res1[2]) %>%
    trimws() %>%
    as.numeric()
  h2_1_se <- sub("p-value", "", res1[3]) %>%
    trimws() %>%
    as.numeric()
  h2_1_pvalue <- res1[4] %>% as.numeric()

  # Ve
  res2 <- strsplit(res[8], ":") %>%
    unlist() %>%
    trimws()
  h2_2 <- sub("SE", "", res2[2]) %>%
    trimws() %>%
    as.numeric()
  h2_2_se <- sub("p-value", "", res2[3]) %>%
    trimws() %>%
    as.numeric()
  h2_2_pvalue <- sub("p-value", "", res2[4]) %>%
    as.numeric()

  # cor
  res3 <- strsplit(res[9], ":") %>%
    unlist() %>%
    trimws()
  cor <- sub("SE", "", res3[2]) %>%
    trimws() %>%
    as.numeric()
  cor_se <- sub("p-value", "", res3[3]) %>%
    trimws() %>%
    as.numeric()
  cor_pvalue <- sub("p-value", "", res3[4]) %>%
    as.numeric()

  save_res <- rbind(
    data.frame(type = "h2_1", value = h2_1, se = h2_1_se, pvalue = h2_1_pvalue),
    data.frame(type = "h2_2", value = h2_2, se = h2_2_se, pvalue = h2_2_pvalue),
    data.frame(type = "cor", value = cor, se = cor_se, pvalue = cor_pvalue)
  )
  return(save_res)
}



splmm <- function(exp, coord, celltype_prop, sel_celltype = NULL, sel_gene = NULL, path_mtg = "./path_mtg", tmpdir = "./tmp", nthread = 1, verbose = 1, remove_tmpdir = TRUE) {
  cat("[", format(Sys.time()), "]", " - Start\n", sep = "")

  if (!dir.exists(tmpdir)) {
    unlink(tmpdir, recursive = TRUE)
    dir.create(tmpdir, recursive = TRUE)
  }

  if (verbose == 2) {
    SYS_PRINT <- FALSE
  } else {
    SYS_PRINT <- TRUE
  }

  dir_current <- getwd()
  setwd(str_glue("{dir_current}/{tmpdir}"))

  ################################################################
  #### 1. preprocessing ##########################################
  ################################################################

  cat("[", format(Sys.time()), "]", " - Calculate spatial kernel / ", sep = "")
  distance_kernel <- cal_spatial_kernel(exp, coord)

  cat("[", format(Sys.time()), "]", " - Calculate celltype kernel (full)\n", sep = "")
  celltype_kernel <- celltype_prop %>% cal_linear_kernel()
  
  cat("[", format(Sys.time()), "]", " - Calculate celltype kernel (sel celltype)\n", sep = "")
  sel_celltype_kernel <- celltype_prop[, colnames(celltype_prop) %in% sel_celltype, drop = FALSE] %>% cal_linear_kernel()
  
  cat("[", format(Sys.time()), "]", " - Calculate spatial kernel (nonsel celltype)\n", sep = "")
  nonsel_celltype_kernel <- celltype_prop[, !colnames(celltype_prop) %in% sel_celltype, drop = FALSE] %>% cal_linear_kernel()

  n_sample <- ncol(exp)
  
  # dummy fam file: ".fam"
  save_fam <- cbind(seq(1, n_sample), seq(1, n_sample), rep(0, n_sample), rep(0, n_sample), -9, -9)
  fwrite(as.data.table(save_fam), "data.fam", sep = "\t", col.names = FALSE)

  # spatial kernel: ".dist"
  fwrite(make_longform(distance_kernel), "kernel.dist", sep = "\t", col.names = FALSE)

  # cell type kernel: ".mat"
  fwrite(make_longform(celltype_kernel), "kernel.mat", sep = "\t", col.names = FALSE)
  fwrite(make_longform(sel_celltype_kernel), "sel_kernel.mat", sep = "\t", col.names = FALSE)
  fwrite(make_longform(nonsel_celltype_kernel), "nonsel_kernel.mat", sep = "\t", col.names = FALSE)

  ###############################################################
  ### 2. Compute kernel matrix Q ################################
  ###############################################################

  cat("[", format(Sys.time()), "]", " - calculate Q matrix (spatial x sel celltype) ", sep = "")
  sel_mat_Q <- cal_Q(distance_kernel, sel_celltype_kernel)
  fwrite(make_longform(sel_mat_Q), "dist_sel.matmat", sep = "\t", col.names = FALSE)

  cat("[", format(Sys.time()), "]", " - calculate Q matrix (spatial x nonsel celltype)", sep = "")
  nonsel_mat_Q <- cal_Q(distance_kernel, nonsel_celltype_kernel)
  fwrite(make_longform(nonsel_mat_Q), "dist_nonsel.matmat", sep = "\t", col.names = FALSE)

  cat("[", format(Sys.time()), "]", " - calculate Q matrix (sel celltype x nonsel celltype)", sep = "")
  sel_nonsel_mat_Q <- cal_Q(sel_celltype_kernel, nonsel_celltype_kernel)
  fwrite(make_longform(sel_nonsel_mat_Q), "sel_nonsel.matmat", sep = "\t", col.names = FALSE)
  
  # 3-0. run greml (null)
  sink("null_greml.matlist")
  cat("kernel.dist", "\n")
  cat("kernel.mat", "\n")
  sink()

  # 3-1. run greml
  sink("greml.matlist")
  cat("kernel.dist", "\n")
  cat("sel_kernel.mat", "\n")
  cat("nonsel_kernel.mat", "\n")
  sink()

  # 3-2. run core greml
  sink("coregreml.matlist")
  cat("kernel.dist", "\n")
  cat("sel_kernel.mat", "\n")
  cat("nonsel_kernel.mat", "\n")
  cat("dist_sel.matmat", "\n")
  cat("dist_nonsel.matmat", "\n")
  cat("sel_nonsel.matmat", "\n")
  sink()
  
  final_res <- pbmcapply::pbmclapply(rownames(exp), function(sel_gene) {
    # expression: ".dat"
    save_exp <- cbind(seq(1, n_sample), seq(1, n_sample), exp[sel_gene, ])
    fwrite(as.data.table(save_exp), str_glue("{sel_gene}_data.dat"), sep = "\t", col.names = FALSE)

    ################################################################
    #### 3. RUN CORE greml #########################################
    ################################################################

    # 3-0. run greml (null)
    cat("[", format(Sys.time()), "]", " - Run GREML\n", sep = "")
    system(str_glue("{dir_current}/{path_mtg} -p data.fam -mg null_greml.matlist -d {sel_gene}_data.dat -mod 1 -thread 1 -out {sel_gene}_null_greml.out"), ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)


    # 3-1. run greml
    cat("[", format(Sys.time()), "]", " - Run GREML\n", sep = "")
    system(str_glue("{dir_current}/{path_mtg} -p data.fam -mg greml.matlist -d {sel_gene}_data.dat -mod 1 -thread 1 -out {sel_gene}_greml.out"), ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)


    # 3-2. run core greml
    cat("[", format(Sys.time()), "]", " - Run CORE GREML\n", sep = "")
    system(str_glue("{dir_current}/{path_mtg} -p data.fam -mg coregreml.matlist -d {sel_gene}_data.dat -mod 1 -thread 1 -out {sel_gene}_coregreml.out"), ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)


    # 3-3. summary variance component
    system(str_glue("grep V {sel_gene}_null_greml.out > {sel_gene}_vc_null_greml"))
    res_null_greml <- fread("{sel_gene}_vc_null_greml")
    res_null_greml[, 1] <- c("v_e", "v_dist", "v_celltype")
    colnames(res_null_greml) <- c("type", "variance", "SE")

    system(str_glue("grep V {sel_gene}_greml.out > {sel_gene}_vc_greml"))
    res_greml <- fread("{sel_gene}_vc_greml")
    res_greml[, 1] <- c("v_e", "v_dist", "v_selCelltype", "v_nonselCelltype")
    colnames(res_greml) <- c("type", "variance", "SE")

    system(str_glue("grep V {sel_gene}_coregreml.out > {sel_gene}_vc_coregreml"))
    res_coregreml <- fread("{sel_gene}_vc_coregreml")
    res_coregreml[, 1] <- c("ve", "v_dist", "v_selCelltype", "v_nonselCelltype", "v_distXselCelltype", "v_distXnonselCelltype", "v_selCelltypeXnonselCelltype")
    colnames(res_coregreml) <- c("type", "variance", "SE")

    # # 3-4. Likelihood-ratio test
    # cat("[", format(Sys.time()), "]", " - Run Likelihood-ratio test\n", sep = "")

    # raw_ll_greml <- system("grep LKH greml.out", intern = TRUE)
    # ll_greml <- sub("LKH", "", raw_ll_greml) %>% trimws() %>% as.numeric()

    # raw_ll_coregreml <- system("grep LKH coregreml.out", intern = TRUE)
    # ll_coregreml <- sub("LKH", "", raw_ll_coregreml) %>% trimws() %>% as.numeric()

    # diff_ll <- (ll_greml - ll_coregreml)
    # res_lrt <- pchisq(-2 * diff_ll, df = 1, lower.tail = FALSE)

    # res_ll <- data.frame(ll_greml = ll_greml, ll_coregreml = ll_coregreml, lrt = diff_ll, pvalue = res_lrt)


    # 4-4. cal heritability (h2)
    cat("[", format(Sys.time()), "]", " - Calculate heritability\n", sep = "")
    if (file.exists(str_glue("{sel_gene}_coregreml.do"))) {
      file.remove(str_glue("{sel_gene}_coregreml.do"))
    }
    system(str_glue("grep -vwE '(LKH|h2)' {sel_gene}_coregreml.out > {sel_gene}_coregreml.out2"))

    sink(str_glue("{sel_gene}_coregreml.do"))
    cat(str_glue("{sel_gene}_coregreml.out2", "\n")) # line 1: specify the file with parameter estimates
    cat(7, "\n") # line 2: tot. # of variance & covariance components in the file
    cat("R 2 1 3 4", "\n") # #2 / (#1+#2+#3+#4)
    cat("R 3 1 2 4", "\n") # #3 / (#1+#2+#3+#4)
    cat("R 4 1 2 3", "\n") # #4 / (#1+#2+#3+#4)
    cat("C 5 2 3", "\n")   # #5 / sqrt(#2 * #3)
    cat("C 6 4 2", "\n")   # #6 / sqrt(#4 * #2)
    cat("C 7 4 3", "\n")   # #7 / sqrt(#4 * #3)
    sink()

    system(str_glue("{dir_current}/{path_mtg} -delta2 {sel_gene}_coregreml.do > res"))
    raw_res <- fread("res", skip = 6, fill = TRUE) %>% as.data.frame()
    
    raw_res_h2 <- raw_res[grep("Ratio", raw_res[, 1]), c(2, 4, 6)]
    colnames(raw_res_h2) <- c("h2", "se", "pvalue")
    res_h2 <- data.frame(
      type = c("dist", "selCelltype", "nonselCelltype"),
      h2 = as.numeric(raw_res_h2$h2),
      se = as.numeric(raw_res_h2$se),
      pvalue = as.numeric(raw_res_h2$pvalue)
    )
    
    raw_res_cor <- raw_res[grep("Cor", raw_res[, 1]), c(3, 5, 7)]
    colnames(raw_res_cor) <- c("cor", "se", "pvalue")
    res_cor <- data.frame(
      type = c("distXselCelltype", "distXnonselCelltype", "selCelltypeXnonselCelltype"),
      cor = as.numeric(raw_res_cor$cor),
      se = as.numeric(raw_res_cor$se),
      pvalue = as.numeric(raw_res_cor$pvalue)
    )
      
    cat("[", format(Sys.time()), "]", " - End\n\n", sep = "")
    setwd(dir_current) 
    
    list(
      nullgreml = res_null_greml,
      greml = res_greml,
      coregrem = res_coregreml,
      # lrt = res_ll,
      h2 = res_h2,
      cor = res_cor
    )
  }, mc.cores = nthread)
  
  if (remove_tmpdir) {
    unlink(str_glue("{tmpdir}"), recursive = TRUE)
  }

  return(final_res)
}
