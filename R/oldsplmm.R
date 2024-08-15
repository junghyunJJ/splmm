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
    cat(paste("## cal bandwidth: ", bandwidth, "\n"), sep = "")
  } else {
    bandwidth <- bandwidth.set.by.user
    cat(paste("## select bandwidth by user: ", bandwidth, "\n"), sep = "")
  }

  # scale coordinate data
  location_normalized <- scale(coord)

  if (sparseKernel == FALSE) {
    cat(paste("## Calculating kernel matrix\n"))
    kernelmat <- SpatialPCA::kernel_build(kerneltype = kerneltype, location = location_normalized, bandwidth = bandwidth)
  } else if (sparseKernel == TRUE) {
    cat(paste("## Calculating sparse kernel matrix\n"))
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
    cat("[", format(Sys.time()), "]", " - Generate near PD: distance kernel\n", sep = "")
    K <- as.matrix(Matrix::nearPD(K)$mat)
  }
  K_chol <- t(chol(K))

  if (!matrixcalc::is.positive.definite(A)) {
    cat("[", format(Sys.time()), "]", " - Generate near PD: celltype kernel\n", sep = "")
    A <- as.matrix(Matrix::nearPD(A)$mat)
  }
  A_chol <- t(chol(A))

  Q <- (A_chol %*% t(K_chol)) + t(A_chol %*% t(K_chol))

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



splmm <- function(exp, distance_kernel, celltype_kernel, path_mtg = "./path_mtg", tmpdir = "./tmp", nthread = 1, verbose = 1, remove_tmpdir = TRUE) {
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

  n_sample <- length(exp)

  # expression: ".dat"
  save_exp <- cbind(seq(1, n_sample), seq(1, n_sample), exp)
  fwrite(as.data.table(save_exp), "data.dat", sep = "\t", col.names = FALSE)

  # dummy fam file: ".fam"
  save_fam <- cbind(seq(1, n_sample), seq(1, n_sample), rep(0, n_sample), rep(0, n_sample), -9, -9)
  fwrite(as.data.table(save_fam), "data.fam", sep = "\t", col.names = FALSE)

  # spatial kernel: ".grm"
  fwrite(make_longform(distance_kernel), "data.grm", sep = "\t", col.names = FALSE)

  # cell type kernel: ".bmat"
  fwrite(make_longform(celltype_kernel), "data.bmat", sep = "\t", col.names = FALSE)


  ###############################################################
  ### 2. Compute kernel matrix Q ################################
  ###############################################################

  cat("[", format(Sys.time()), "]", " - calculate Q matrix\n", sep = "")
  mat_Q <- cal_Q(distance_kernel, celltype_kernel)
  fwrite(make_longform(mat_Q), "grm_bmat.chol.matmat2", sep = "\t", col.names = FALSE)


  ################################################################
  #### 3. RUN CORE greml #########################################
  ################################################################

  # 3-1. run greml
  sink("greml.matlist")
  cat("data.grm", "\n")
  cat("data.bmat", "\n")
  sink()

  cat("[", format(Sys.time()), "]", " - Run GREML\n", sep = "")
  system(str_glue("{dir_current}/{path_mtg} -p data.fam -mg greml.matlist -d data.dat -mod 1 -thread {nthread} -out greml.out"), ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)

  # 3-2. run core greml
  sink("coregreml.matlist")
  cat("data.grm", "\n")
  cat("data.bmat", "\n")
  cat("grm_bmat.chol.matmat2", "\n")
  sink()

  cat("[", format(Sys.time()), "]", " - Run CORE GREML\n", sep = "")
  system(str_glue("{dir_current}/{path_mtg} -p data.fam -mg coregreml.matlist -d data.dat -mod 1 -thread {nthread} -out coregreml.out"), ignore.stdout = SYS_PRINT, ignore.stderr = SYS_PRINT)

  # 3-3. summary variance component
  system("grep V greml.out > vc_greml")
  vc_greml <- fread("vc_greml")
  vc_greml[, 1] <- c("ve", "v1", "v2")
  colnames(vc_greml) <- c("type", "variance", "SE")
  res_greml <- cbind(type = "greml", vc_greml)

  system("grep V coregreml.out > vc_coregreml")
  vc_coregreml <- fread("vc_coregreml")
  vc_coregreml[, 1] <- c("ve", "v1", "v2", "v12")
  colnames(vc_coregreml) <- c("type", "variance", "SE")
  res_coregreml <- cbind(type = "coregreml", vc_coregreml)

  # 3-4. Likelihood-ratio test
  cat("[", format(Sys.time()), "]", " - Run Likelihood-ratio test\n", sep = "")

  raw_ll_greml <- system("grep LKH greml.out", intern = TRUE)
  ll_greml <- sub("LKH", "", raw_ll_greml) %>% trimws() %>% as.numeric()

  raw_ll_coregreml <- system("grep LKH coregreml.out", intern = TRUE)
  ll_coregreml <- sub("LKH", "", raw_ll_coregreml) %>% trimws() %>% as.numeric()

  diff_ll <- (ll_greml - ll_coregreml)
  res_lrt <- pchisq(-2 * diff_ll, df = 1, lower.tail = FALSE)

  res_ll <- data.frame(ll_greml = ll_greml, ll_coregreml = ll_coregreml, lrt = diff_ll, pvalue = res_lrt)

  # 4-4. cal heritability (h2)
  cat("[", format(Sys.time()), "]", " - Calculate heritability\n", sep = "")
  if (file.exists("coregreml.do")) {
    file.remove("coregreml.do")
  }
  system("grep -vwE '(LKH|h2)' coregreml.out > coregreml.out2")

  sink("coregreml.do")
  cat("coregreml.out2", "\n") # line 1: specify the file with parameter estimates
  cat(4, "\n") # line 2: tot. # of variance & covariance components in the file
  cat("R 2 1 3 4 4", "\n") # line 3: compute prop. of variance due to genetics
  cat("R 3 1 2 4 4", "\n") # line 4: compute prop. of variance due to environments
  cat("C 4 2 3", "\n") # line 5: compute correlation between g & b
  sink()

  raw_res_h2 <- system(str_glue("{dir_current}/{path_mtg} -delta2 coregreml.do"), intern = TRUE)
  res_h2 <- summry_res(raw_res_h2)

  final_res <- list(
    greml = res_greml,
    coregrem = res_coregreml,
    lrt = res_ll,
    h2 = res_h2
  )

  cat("[", format(Sys.time()), "]", " - End\n\n", sep = "")
  setwd(dir_current)

  if (remove_tmpdir) {
    unlink(str_glue("{tmpdir}"), recursive = TRUE)
  }

  return(final_res)
}
