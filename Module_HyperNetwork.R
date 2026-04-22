rm(list = ls())
load("funclu.RData")

library(igraph)
library(pracma)
library(ADSIHT)
library(Matrix)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(patchwork)
get_df.fit <- function(){
  times_a = times[1:3];times_b = times[4:6];times_c = times[7:9];times_d = times[10:12];
  times_aa = seq(min(times_a),max(times_a),length = 30)
  times_bb = seq(min(times_b),max(times_b),length = 30)
  times_cc = seq(min(times_c),max(times_c),length = 30)
  times_dd = seq(min(times_d),max(times_d),length = 30)
  
  par.mu = result$mu_par
  k = result$cluster_number
  
  mu.fit = list(power_equation(times_aa, par.mu[1:k,1:2]), 
                power_equation(times_bb, par.mu[1:k,3:4]),
                power_equation(times_cc, par.mu[1:k,5:6]),
                power_equation(times_dd, par.mu[1:k,7:8]))
  for (i in 1:4) {
    rownames(mu.fit[[i]]) = paste0('M',1:18)
  }
  return(mu.fit)
}

df.fit = get_df.fit()

module.df = lapply(1:18, function(xi)
  rbind(df.fit[[1]][xi,],df.fit[[2]][xi,],df.fit[[3]][xi,], df.fit[[4]][xi,]))



poly_basis_1d <- function(x, name = NULL, degree = 3) {
  n <- length(x)
  X <- matrix(NA, n, degree)
  for (d in 1:degree) {
    X[, d] <- x^d
  }
  
  base_name <- if (!is.null(name)) name else "x"
  colnames(X) <- paste0(base_name, "_", 1:degree)
  
  return(X)
}






MTODE <- function(Y, times, M = 5, smooth = "bs",effect_thr = 1e-5,name = 'NULL') {
  if( is.null(colnames(Y))){
    colnames(Y) = paste0("x",1:ncol(Y))
  }
  
  n <- length(times)
  p <- ncol(Y)
  
  original_min <- min(times)
  original_range <- max(times) - min(times)
  # 归一化时间
  times_norm <- (times - min(times)) / (max(times) - min(times))
  times_new = seq(min(times_norm),max(times_norm),length = 50)
  times_restored <- times_new * original_range + original_min
  
  # 1. 预平滑 (Smoothing Splines)
  if (smooth == "bs") {
    fit_spline <- lapply(1:p, function(xi) smooth.spline(times_norm, as.numeric(Y[, xi])))
    x_smooth <- sapply(1:p, function(xi) predict(fit_spline[[xi]], x = times_norm)$y)
    colnames(x_smooth) <- colnames(Y)
    
    x_smooth2 <- sapply(1:p, function(xi) predict(fit_spline[[xi]], x = times_new)$y)
    colnames(x_smooth2) <- colnames(Y)
    
  } else if(smooth == "power_equation")
  {
    
    times_safe <- ifelse(times == 0, 1e-6, times)
    fit_results <- data.frame(
      Variable  = colnames(Y),   # 如果 Y 没有列名，会自动生成
      a_final   = numeric(p),
      b_final   = numeric(p),
      is_nls    = logical(p), # 记录 nls 是否成功收敛
      stringsAsFactors = FALSE
    )
    
    x_smooth <- matrix(NA, nrow = length(times_norm), ncol = ncol(Y))
    colnames(x_smooth) <- colnames(Y)
    x_smooth2 <- matrix(NA, nrow = length(times_new), ncol = ncol(Y))
    colnames(x_smooth2) <- colnames(Y)
    
    for (i in 1:p) {
      
      y <- Y[, i]
      y_safe <- ifelse(y == 0, 1e-6, y)
      
      fit_init <- try(lm(log(y_safe) ~ log(times_safe)), silent = TRUE)
      
      if (inherits(fit_init, "try-error")) {
        # 如果极度异常，全部填 NA 并跳过本次循环
        fit_results$a_final[i] <- NA
        fit_results$b_final[i] <- NA
        fit_results$is_nls[i] <- FALSE
        next
      }
      
      a_init <- as.numeric(exp(coef(fit_init)[1]))
      b_init <- as.numeric(coef(fit_init)[2])
      
      fit_nls <- tryCatch({
        nls(y ~ a * times^b, 
            start = list(a = a_init, b = b_init),
            control = nls.control(maxiter = 200, warnOnly = TRUE)) # warnOnly允许未完美收敛时也返回结果
      }, error = function(e) {
        # 如果 nls 彻底报错失败，返回 NULL
        return(NULL) 
      })
      
      # 第三步：结果记录
      if (!is.null(fit_nls)) {
        # nls 成功：记录精准参数
        fit_results$a_final[i] <- coef(fit_nls)["a"]
        fit_results$b_final[i] <- coef(fit_nls)["b"]
        fit_results$is_nls[i]  <- TRUE
        a_fit <- coef(fit_nls)["a"]
        b_fit <- coef(fit_nls)["b"]
        if (!is.na(a_fit) && !is.na(b_fit)) {
          # 将原始 times 代入模型：y = a * x^b
          x_smooth[, i] <- a_fit * (times_norm ^ b_fit)
          x_smooth2[, i] <- a_fit * (times_new ^ b_fit)
          
        } else {
          x_smooth[, i] <- NA # 或者换成 Y[, i] 保留原值
          x_smooth2[, i] <- NA # 或者换成 Y[, i] 保留原值
        }
        
      } else {
        # nls 失败：退而求其次，记录第一步 lm 获取的近似参数
        fit_results$a_final[i] <- a_init
        fit_results$b_final[i] <- b_init
        fit_results$is_nls[i]  <- FALSE
      }
      
    }
  }
  
  
  # 2. 构建基函数并积分
  # 构造多项式基
  phi_list <- lapply(1:p, function(xi) {
    poly_basis_1d(x_smooth[, xi], name = colnames(Y)[xi], degree = M)
  })
  phi <- Reduce(cbind, phi_list)
  
  # 对基进行积分 (Trapezoidal integration)
  phi_int <- apply(phi, 2, function(col) cumtrapz(times_norm, col))
  
  # --- 定义组构造函数 ---
  construct_group <- function(j, phi_int_data) {
    s0 <- p * (j - 1) + 1
    group <- rep(s0:(s0 + p - 1), each = M)
    group0 <- min(group)
    group <- sort(c(group0, group))
    group_levels <- unique(group)
    
    phi1_int_local <- cbind(times_norm, phi_int_data)
    group_idx_list <- split(seq_len(ncol(phi1_int_local)), group)
    
    qr_list <- lapply(group_levels, function(g) {
      cols <- group_idx_list[[as.character(g)]]
      Xg <- phi1_int_local[, cols, drop = FALSE]
      mXg <- colMeans(Xg)
      Xg_centered <- sweep(Xg, 2, mXg, FUN = "-")
      
      # 若秩不足则加微小噪声
      if (qr(Xg_centered)$rank < ncol(Xg_centered)) {
        Xg_centered <- Xg_centered + matrix(rnorm(length(Xg_centered), sd = 1e-12),
                                            nrow = nrow(Xg_centered), ncol = ncol(Xg_centered))
      }
      qr_g <- qr(Xg_centered)
      list(Q = qr.Q(qr_g, complete = FALSE), R = qr.R(qr_g), mXg = mXg, src_cols = cols)
    })
    
    Q_all <- do.call(cbind, lapply(qr_list, `[[`, "Q"))
    
    list(
      group = group,
      group_idx_list = group_idx_list,
      phi_int = phi_int_data,
      Q_all = Q_all,
      R_list = lapply(qr_list, `[[`, "R"),
      mXg_list = lapply(qr_list, `[[`, "mXg"),
      src_cols_list = lapply(qr_list, `[[`, "src_cols")
    )
  }
  
  # 3. 构造大设计矩阵 (Block Diagonal Q)
  Q_info <- lapply(1:p, construct_group, phi_int_data = phi_int)
  Q_blocks <- lapply(Q_info, `[[`, "Q_all")
  X_design <- as.matrix(bdiag(Q_blocks))
  
  # 4. 数据标准化与回归
  y_scaled <- scale(Y)
  Y_all <- matrix(as.vector(y_scaled), ncol = 1)
  
  groups_vec <- unlist(lapply(Q_info, `[[`, "group"))
  
  # ADSIHT 回归
  # fit <- ADSIHT(X_design, Y_all, group = groups_vec)
  
  fit <- ADSIHT(X_design, Y_all, group = groups_vec)
  
  
  best_idx <- which.min(fit$ic)
  beta_std_Q_vec <- fit$beta[, best_idx]
  
  #fit <- cv.sparsegl(X_design, Y_all, group = groups_vec) #sparsegl
  #fit <- cv.grpreg(X_design, Y_all, group = groups_vec)   #grpreg
  #beta_std_Q_vec <- as.vector(coef(fit, s = fit$lambda.1se)[-1])
  
  # 5. 恢复原始 Beta
  recover_beta <- function(beta_std_Q_vec, Q_info, y_s) {
    cols_per_j <- sapply(Q_info, function(z) ncol(z$Q_all))
    cum_cols <- c(0, cumsum(cols_per_j))
    sY <- attr(y_s, "scaled:scale")
    mY <- attr(y_s, "scaled:center")
    
    res_list <- lapply(seq_along(Q_info), function(j) {
      idx_range <- (cum_cols[j] + 1):cum_cols[j + 1]
      bQj_all <- beta_std_Q_vec[idx_range]
      
      info <- Q_info[[j]]
      rcols_vec <- sapply(info$R_list, ncol)
      bQg_list <- split(bQj_all, rep(seq_along(rcols_vec), rcols_vec))
      
      beta_Xj_full <- numeric(ncol(info$phi_int))
      mX_full <- numeric(ncol(info$phi_int))
      
      Map(function(Rg, bQg, src_cols, mXg) {
        # 注意: backsolve 结果写入父环境变量
        beta_Xg <- backsolve(Rg, bQg)
        beta_Xj_full[src_cols] <<- beta_Xg
        mX_full[src_cols] <<- mXg
      }, info$R_list, bQg_list, info$src_cols_list, info$mXg_list)
      
      beta_orig_j <- sY[j] * beta_Xj_full
      intercept_orig_j <- mY[j] - sY[j] * sum(mX_full * beta_Xj_full)
      list(beta = beta_orig_j, intercept = intercept_orig_j)
    })
    
    list(
      B_est = lapply(res_list, `[[`, "beta"),
      intercept = sapply(res_list, `[[`, "intercept")
    )
  }
  
  beta_recovered <- recover_beta(beta_std_Q_vec, Q_info, y_scaled)
  beta_orig <- Reduce(cbind, beta_recovered$B_est)
  
  beta_hat <- beta_orig[-1, ]
  beta_all <- rbind(beta_recovered$intercept, beta_orig)
  
  
  
  # 6. 计算估计值与绘图
  #new_phi_all =  cbind(1, as.matrix(bdiag(cbind(times_norm, phi_int))))
  phi_all <- cbind(1, as.matrix(bdiag(cbind(times_norm, phi_int))))
  y_est <- phi_all %*% beta_all
  
  y_est_df <- melt(data.frame(times = times, y_est), id.vars = 'times')
  y_obs_df <- melt(data.frame(times = times, Y), id.vars = 'times')
  
  pp0 <- ggplot() + 
    geom_point(data = y_obs_df, aes(times, value, group = variable)) +
    geom_line(data = y_est_df, aes(times, value, group = variable), color = 'red') +
    facet_wrap(~variable) + 
    theme_bw()
  
  
  
  
  # 找出固定为0的参数索引
  par_init <- c(beta_orig)        # 初始参数向量
  par_mat <- matrix(par_init, ncol = ncol(Y))
  fixed_idx <- which(par_init == 0)
  free_idx  <- setdiff(seq_along(par_init), fixed_idx)
  
  # ridge目标函数
  # 1. 准备数据
  lambda <- 1e-1
  Y_c <- as.matrix(sweep(Y, 2, beta_recovered$intercept, "-"))
  X <- as.matrix(cbind(times_norm, phi_int))
  
  pp <- ncol(X)     # 特征总数 (对应 par_mat 的行数)
  m <- ncol(Y_c)   # 响应变量数 (对应 par_mat 的列数)
  
  # 2. 初始化最终的参数矩阵 (全部置 0)
  par_mat <- matrix(0, nrow = pp, ncol = m)
  
  # 3. 按列 (Y 的每一个变量) 进行解析求解
  for (j in 1:m) {
    
    start_idx <- (j - 1) * pp + 1
    end_idx   <- j * pp
    
    global_free_for_j <- free_idx[free_idx >= start_idx & free_idx <= end_idx]
    local_free <- global_free_for_j - start_idx + 1
    
    if (length(local_free) > 0) {
      
      X_sub <- X[, local_free, drop = FALSE]
      y_sub <- Y_c[, j]
      
      # 构建单位矩阵 I
      I <- diag(length(local_free))
      
      # 构建 A 矩阵和 b 向量 (Ax = b 形式)
      A <- t(X_sub) %*% X_sub + lambda * I
      b <- t(X_sub) %*% y_sub
      
      # ---------------------------------------------------------
      # 核心修复：双重防弹求逆机制
      # ---------------------------------------------------------
      beta_ridge <- tryCatch({
        solve(A, b)
        
      }, error = function(e) {
        MASS::ginv(A) %*% b
      })
      
      # 填回参数
      par_mat[local_free, j] <- as.numeric(beta_ridge)
    }
  }
  
  # (可选) 将参数矩阵展平，与之前 optim 的格式保持一致
  opt_par_full <- as.numeric(par_mat)
  final_free_par <- opt_par_full[free_idx]
  
  
  
  # 得到最终参数（原本为0的保持0）
  par_final <- par_init
  par_final[free_idx] <- final_free_par
  par_final_mat <- matrix(par_final, ncol = ncol(Y))
  beta_all = rbind(beta_recovered$intercept,par_final_mat)
  
  
  phi_list2 <- lapply(1:p, function(xi) {
    poly_basis_1d(x_smooth2[, xi], name = colnames(Y)[xi], degree = M)
  })
  phi2 <- Reduce(cbind, phi_list2)
  
  # 对基进行积分 (Trapezoidal integration)
  phi_int2 <- apply(phi2, 2, function(col) cumtrapz(times_new, col))
  phi_all2 <- cbind(1, as.matrix(bdiag(cbind(times_new, phi_int2))))
  
  
  y_est2 <- phi_all2 %*% beta_all
  
  y_est_df2 <- melt(data.frame(times = times_restored, y_est2), id.vars = 'times')
  y_obs_df <- melt(data.frame(times = times, Y), id.vars = 'times')
  
  pp1 <- ggplot() + 
    geom_point(data = y_obs_df, aes(times, value, group = variable)) +
    geom_line(data = y_est_df2, aes(times, value, group = variable), color = 'red') +
    facet_wrap(~variable) + 
    theme_bw()
  
  
  
  
  
  #绘制残差
  residuals_mat <- Y - phi_all %*% beta_all
  colnames(residuals_mat) <- colnames(Y) # 继承变量名
  
  # 将残差矩阵转换为长数据格式以适配 ggplot
  resid_df <- melt(data.frame(times = times_norm, residuals_mat), id.vars = 'times')
  p_resid_time <- ggplot(resid_df, aes(x = times, y = value, color = variable)) +
    # 画出残差散点
    geom_point(size = 1.5, alpha = 0.8) +
    # 画出连接到 y=0 的垂直线，形成“棒棒糖”视觉效果
    geom_segment(aes(x = times, xend = times, y = 0, yend = value), alpha = 0.5) +
    # 添加一条基准的零线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_wrap(~variable, scales = "free_y") + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none") + # 分面已经自带变量信息，可隐藏图例
    labs(title = "Residuals vs. Time",
         x = "Times", 
         y = "Observed - Estimated")
  
  
  
  
  
  
  # 计算组贡献 (f_group_est)
  f_group <- function(jj) {
    
    betam <- beta_all[-c(1:2), jj]
    est <- matrix(0, nrow = length(times_norm), ncol = p)
    
    group_idx_list <- split(seq_len(ncol(phi)), rep(1:p, each = M))
    
    for (g in seq_len(p)) {
      est[, g] <- phi_int[, group_idx_list[[g]], drop = FALSE] %*% betam[group_idx_list[[g]]]
    }
    
    # 加上截距项和时间线性项
    est[, jj] <- est[, jj] + beta_all[2, jj] * times_norm +  beta_all[1, jj]
    
    colnames(est) <- colnames(Y)
    return(est)
  }
  
  f_group_est <- lapply(1:p, f_group)
  
  
  
  get_active_HOI <- function(){
    get_possible_edge <- function(j) {
      edge <- NULL 
      names <- c("H1", "H2", "S1", "S2")
      for (s in (1:4)[-j]) {
        # Add drop = FALSE to prevent demotion of a single column and loss of column names
        tmp_mat <- f_group_est[[s]][, -c(j, s), drop = FALSE] 
        
        # Modify the column names
        colnames(tmp_mat) <- paste0(names[s], "<-", colnames(tmp_mat))
        
        if (is.null(edge)) {
          edge <- tmp_mat
        } else {
          edge <- cbind(edge, tmp_mat) 
        }
      }
      col_sums_abs <- colSums(abs(edge))
      edge <- edge[, col_sums_abs != 0, drop = FALSE]
      return(edge)
    }
    
    edges <- lapply(1:4, function(xi) get_possible_edge(xi))
    
    edges <- lapply(edges, function(df) {
      # 【防错1】：如果传入的df已经是0列，直接返回
      if (ncol(df) == 0) return(df)
      
      col_means <- colMeans(abs(df))
      keep_cols <- names(col_means)[col_means >= effect_thr*10]
      df_filtered <- df[, keep_cols, drop = FALSE]
      
      return(df_filtered)
    })
    
    num_blocks <- length(edges) 
    n_times <- length(times_norm) # 获取时间点长度，用于构建占位空矩阵
    
    # 1. Iterate through the 4 elements
    phi_high_int_list <- lapply(1:num_blocks, function(i) {
      
      # 【防错2】：如果当前组没有边，返回一个 n_times 行 0 列的矩阵
      # 这极其重要！它能保证后续 bdiag() 拼接时，这一块依然保留行空间，不会导致总行数错乱
      if (ncol(edges[[i]]) == 0) {
        return(matrix(0, nrow = n_times, ncol = 0))
      }
      
      # 【防错3】：使用 seq_len(ncol) 替代 1:ncol，防止出现 1:0 的致命 bug
      phi_list_high <- lapply(seq_len(ncol(edges[[i]])), function(xi) {
        poly_basis_1d(edges[[i]][, xi], name = colnames(edges[[i]])[xi], degree = M)
      })
      
      # Merge by column
      phi_high <- Reduce(cbind, phi_list_high)
      
      # Integrate the basis
      phi_high_int <- apply(phi_high, 2, function(col) cumtrapz(times_norm, col))
      
      return(phi_high_int)
    })
    
    # 2. Assemble the block diagonal matrix
    phi_high_int_bdiag <- as.matrix(bdiag(phi_high_int_list))
    
    # 【防错4】：统计全局所有边数。如果为0，直接跳过 ADSIHT，防止报错
    total_edges <- sum(sapply(edges, ncol))
    error <- as.vector(residuals_mat)
    
    if (total_edges > 0) {
      # 只有存在有效边时，才构建 group 并进行拟合
      high_group <- rep(1:total_edges, each = M)
      
      fit_high <- ADSIHT(phi_high_int_bdiag, error, group = high_group)
      beta_high <- fit_high$beta[, which.min(fit_high$ic)]
      
      free_idx_high <- which(beta_high != 0)
      beta_high_final <- beta_high
      
    } else {
      # 如果四组全被过滤空了，给变量赋初始空值，安全通过
      cat("--- 警告：未发现任何有效边，跳过 ADSIHT 拟合 ---\n")
      free_idx_high <- integer(0)
      beta_high_final <- numeric(0)
    }
    
    
    if (length(free_idx_high) > 0) {
      
      # Extract the valid feature subset (drop = FALSE is mandatory to prevent demotion error for a single column)
      X_sub_high <- phi_high_int_bdiag[, free_idx_high, drop = FALSE]
      y_sub_high <- error
      
      lambda_high <- 1e-5
      
      # Construct identity matrix I
      I_high <- diag(length(free_idx_high))
      
      # Construct matrix A and vector b (Ax = b format)
      A_high <- t(X_sub_high) %*% X_sub_high + lambda_high * I_high
      b_high <- t(X_sub_high) %*% y_sub_high
      
      beta_ridge_high <- tryCatch({
        # Attempt 1: Standard linear system solver
        solve(A_high, b_high)
      }, error = function(e) {
        # Attempt 2: If extreme collinearity of polynomial bases causes singularity, use SVD generalized inverse
        MASS::ginv(A_high) %*% b_high
      })
      
      # 4. Fill the Ridge re-calibrated coefficients back into the original vector
      beta_high_final[free_idx_high] <- as.numeric(beta_ridge_high)
      
      cat("--- High-order network Ridge refitting successful! ---\n")
      cat("Activated and optimized", length(free_idx_high), "high-order polynomial parameters.\n")
      
    } else {
      cat("--- No significant high-order interaction effects found, skipping Ridge fitting. ---\n")
    }
    
    
    
    
    beta_high <- beta_high_final
    
    # Recalculate fitted errors using optimized parameters
    error_hat_vec <- as.numeric(phi_high_int_bdiag %*% beta_high)
    residuals_hat_mat <- matrix(error_hat_vec, 
                                nrow = nrow(residuals_mat), 
                                ncol = ncol(residuals_mat))
    
    # Inherit the column names of the original residual matrix (i.e., names of the 4 nodes)
    colnames(residuals_hat_mat) <- colnames(residuals_mat)
    
    df_error_obs <- melt(data.frame(times = times_norm, residuals_mat), 
                         id.vars = "times", variable.name = "Node", value.name = "Observed_Error")
    
    df_error_fit <- melt(data.frame(times = times_norm, residuals_hat_mat), 
                         id.vars = "times", variable.name = "Node", value.name = "Fitted_Error")
    
    num_nodes <- ncol(residuals_mat)
    node_names <- colnames(residuals_mat)
    
    df_totals_list <- list()
    df_effects_list <- list()
    
    # Used to track the current index in the global beta_high vector
    current_idx <- 1 
    
    for (i in 1:num_nodes) {
      
      node_name <- node_names[i]
      edge_mat <- edges[[i]]
      num_edges <- ncol(edge_mat)
      
      # If there are no column names, automatically generate fallback names
      edge_names <- colnames(edge_mat)
      if (is.null(edge_names)) edge_names <- paste0("Edge_", 1:num_edges)
      
      # Total number of parameters corresponding to the i-th node
      num_params <- num_edges * M
      
      # [BUG FIX]: Safe extraction to avoid "1:0" reverse sequence bug if num_params is 0
      if (num_params > 0) {
        # Accurately extract the parameter block belonging to the current node from the global beta_high
        beta_i <- beta_high[current_idx : (current_idx + num_params - 1)]
        phi_i <- phi_high_int_list[[i]]
        
        # Calculate the total high-order fitted effect for this node (thick red line)
        total_fit <- phi_i %*% beta_i
      } else {
        beta_i <- numeric(0)
        total_fit <- rep(0, length(times_norm))
      }
      
      df_totals_list[[i]] <- data.frame(
        times = as.numeric(times_norm),
        Node = node_name,
        Observed = as.numeric(residuals_mat[, i]),
        Fitted = as.numeric(total_fit)
      )
      
      # Further decomposition: iterate through all candidate edges under this node
      effects_mat <- matrix(NA, nrow = length(times_norm), ncol = 0)
      valid_edge_names <- c()
      
      if (num_edges > 0) {
        for (k in 1:num_edges) {
          # Extract the M parameters of the k-th edge
          idx_k <- ((k - 1) * M + 1) : (k * M)
          beta_k <- beta_i[idx_k]
          
          # [Core Logic]: If these M parameters are not all 0, it means this high-order edge was selected by the algorithm!
          if (sum(abs(beta_k)) > 0) {
            # Calculate the independent high-order effect contributed by this edge
            effect_k <- phi_i[, idx_k, drop = FALSE] %*% beta_k
            effects_mat <- cbind(effects_mat, as.numeric(effect_k))
            valid_edge_names <- c(valid_edge_names, edge_names[k])
          }
        }
      }
      
      # Convert the effects of non-zero edges into long data format
      if (ncol(effects_mat) > 0) {
        colnames(effects_mat) <- valid_edge_names
        df_eff <- melt(data.frame(times = as.numeric(times_norm), effects_mat), 
                       id.vars = "times", variable.name = "Edge", value.name = "Effect")
        df_eff$Node <- node_name
        df_effects_list[[i]] <- df_eff
      }
      
      # Update the index and move to the next node
      current_idx <- current_idx + num_params
    }
    
    
    
    
    
    # Merge data from all nodes
    df_totals <- do.call(rbind, df_totals_list)
    if(length(df_effects_list) > 0) {
      df_effects <- do.call(rbind, df_effects_list)
    } else {
      df_effects <- data.frame(times=numeric(0), Edge=character(0), Effect=numeric(0), Node=character(0))
    }
    
    # ==========================================
    # 2. Extract the endpoint of each high-order effect line for smart end labeling
    # ==========================================
    if(nrow(df_effects) > 0) {
      max_t <- max(df_effects$times, na.rm = TRUE)
      df_labels <- df_effects[df_effects$times == max_t, ]
    } else {
      df_labels <- data.frame()
    }
    
    # ==========================================
    # 3. Plot the high-order effect decomposition grid for the 4 subplots
    # ==========================================
    p_hoi_decomp <- ggplot() +
      
      # 1. Background baseline
      geom_hline(yintercept = 0, color = "darkgray", linewidth = 0.8) +
      
      # 2. Plot the independent effect of each non-zero edge (colored dashed lines)
      geom_line(data = df_effects, aes(x = times, y = Effect, color = Edge), 
                linewidth = 0.8, linetype = "dashed", alpha = 0.8) +
      
      # 3. Plot the overall high-order fitted error (thick solid red line)
      geom_line(data = df_totals, aes(x = times, y = Fitted), 
                color = "firebrick", linewidth = 1.2) +
      
      # 4. Original residual data points (black semi-transparent points)
      geom_point(data = df_totals, aes(x = times, y = Observed), 
                 color = "black", size = 1.5, alpha = 0.5) +
      
      # 5. [Direct End Labeling] Avoid looking at the legend for colors
      geom_text_repel(data = df_labels, 
                      aes(x = times, y = Effect, label = Edge, color = Edge),
                      nudge_x = 0.05 * diff(range(as.numeric(times_norm))),
                      direction = "y", hjust = 0, segment.size = 0.3,
                      fontface = "bold", size = 3.5, show.legend = FALSE) +
      
      # 6. Facet display by Node
      facet_wrap(~ Node, scales = "free_y", ncol = 2) +
      
      # Leave space on the right for labels
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) + 
      
      # Overall beautification
      theme_bw(base_size = 14) +
      theme(
        legend.position = "none", # Since we have end labels, hide the redundant legend
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "High-Order Network Decomposition of Residuals (HOI Additive Decomposition)",
        subtitle = "Solid red line (Total HOI effect) = Sum of colored dashed lines (HOI contribution of specific edges)",
        x = "Normalized Time",
        y = "High-order Dynamic Effect (Effect Size)"
      )
    
    # Render the plot
    print(p_hoi_decomp)
    
    
    
    
    
    
    # ==========================================
    # 4. 绘制 plot_hoi_2：指定分面顺序与极简透明标签
    # ==========================================
    require(ggrepel)
    require(reshape2)
    require(ggplot2)
    
    # 1. 整理观测数据 (散点：原始 Y 数据)
    df_obs <- melt(data.frame(times = times, Y), id.vars = "times", 
                   variable.name = "Node", value.name = "Observed")
    
    # 2. 整理 Independent, Pairwise(Dependent) 和 总拟合 数据
    df_pairwise_list <- list()
    df_ind_list <- list()
    df_total_fit_list <- list()
    
    for (i in 1:p) {
      node_name <- colnames(Y)[i]
      pw_mat <- f_group_est[[i]]
      
      # a. 提取独立效应 (Independent effect)
      ind_effect <- pw_mat[, i]
      df_ind_list[[i]] <- data.frame(
        times = times,
        Node = node_name,
        Effect = ind_effect,
        Edge = "independent"
      )
      
      # b. 提取依赖效应 (Dependent/Pairwise effects)
      dep_mat <- pw_mat[, -i, drop = FALSE]
      if (ncol(dep_mat) > 0) {
        # 1. 先赋好列名 (边名称)
        colnames(dep_mat) <- paste0(colnames(Y)[-i], "->", node_name)
        
        # 2. 【核心修改】：计算每条边在所有时间点上的绝对值之和，过滤掉全为 0 的列
        # (使用 > 0 确保精确剔除没有任何动态效应的冗余边)
        keep_edges <- colSums(abs(dep_mat)) > 0 
        dep_mat <- dep_mat[, keep_edges, drop = FALSE]
        
        # 3. 只有在过滤后依然存在有效边的情况下，才进行宽表转长表的操作
        if (ncol(dep_mat) > 0) {
          df_pw <- melt(data.frame(times = times, dep_mat), id.vars = "times",
                        variable.name = "Edge", value.name = "Effect")
          df_pw$Node <- node_name
          df_pairwise_list[[i]] <- df_pw
        }
      }
      
      # c. 计算总和
      pw_total <- rowSums(pw_mat)
      hoi_subset <- df_totals[df_totals$Node == node_name, ]
      hoi_total <- if (nrow(hoi_subset) > 0) hoi_subset$Fitted else rep(0, length(times))
      
      # 总拟合
      total_fit <- pw_total + hoi_total
      
      df_total_fit_list[[i]] <- data.frame(
        times = times,
        Node = node_name,
        Total_Fit = total_fit
      )
    }
    
    
    
    
    df_pairwise <- do.call(rbind, df_pairwise_list)
    df_ind <- do.call(rbind, df_ind_list)
    df_total_fit <- do.call(rbind, df_total_fit_list)
    
    # 3. 整理 HOI 效应
    if (exists("df_effects") && nrow(df_effects) > 0) {
      df_hoi <- df_effects
      df_hoi$times <- df_hoi$times * original_range + original_min
    } else {
      df_hoi <- data.frame(times = numeric(0), Edge = character(0), Effect = numeric(0), Node = character(0))
    }
    
    # ==========================================
    # 【新增核心逻辑】：强制锁定分面顺序
    # ==========================================
    target_order <- c("S1", "S2", "H1", "H2")
    
    # 将所有绘图数据框的 Node 转换为因子，并指定 levels
    df_obs$Node <- factor(df_obs$Node, levels = target_order)
    df_pairwise$Node <- factor(df_pairwise$Node, levels = target_order)
    df_ind$Node <- factor(df_ind$Node, levels = target_order)
    df_total_fit$Node <- factor(df_total_fit$Node, levels = target_order)
    
    if (nrow(df_hoi) > 0) {
      df_hoi$Node <- factor(df_hoi$Node, levels = target_order)
    }
    
    # ==========================================
    # 提取末端标签数据
    # ==========================================
    max_t_orig <- max(times, na.rm = TRUE)
    
    df_pw_labels <- df_pairwise[df_pairwise$times == max_t_orig, ]
    df_ind_labels <- df_ind[df_ind$times == max_t_orig, ]
    
    if (nrow(df_hoi) > 0) {
      max_t_hoi <- max(df_hoi$times, na.rm = TRUE)
      df_hoi_labels <- df_hoi[abs(df_hoi$times - max_t_hoi) < 1e-6, ] 
    } else {
      df_hoi_labels <- data.frame(times = numeric(0), Edge = character(0), Effect = numeric(0), Node = character(0))
    }
    
    
    df_pw_labels$Edge <- substr(df_pw_labels$Edge, 1, 2)
    df_hoi_labels$Edge <- gsub("..", "<-", df_hoi_labels$Edge, fixed = TRUE)
    
    
    
    # ==========================================
    # 绘制最终的 plot_hoi_2
    # ==========================================
    bg_df <- data.frame(
      Node = factor(target_order, levels = target_order),
      fill  = c("#FDECEA", "#FDECEA", "#E8F4FC", "#E8F4FC")
    )
    
    plot_hoi_2 <- ggplot() +
      # 背景色：S1/S2 浅红，H1/H2 浅蓝
      geom_rect(data = bg_df,
                aes(fill = fill),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                alpha = 1, inherit.aes = FALSE) +
      scale_fill_identity() +
      # f. 原始观测数据
      geom_point(data = df_obs, aes(x = times, y = Observed),
                 color = "#53CBF3", size = 1.5, alpha = 0.4) +
      
      # a. 零线基准
      geom_hline(yintercept = 0, linetype = 2, color = "darkgray") +
      
      # b. Pairwise 效应 (Dependent) - 绿色 #67AE6E
      geom_line(data = df_pairwise, aes(x = times, y = Effect, group = Edge), 
                color = "#67AE6E", linewidth = 0.8, alpha = 0.8) +
      
      # c. Independent 效应 (自身) - 红色 #DC2525
      geom_line(data = df_ind, aes(x = times, y = Effect, group = Edge), 
                color = "#DC2525", linewidth = 1) +
      
      # d. 互作效应 HOI - 紫色
      geom_line(data = df_hoi, aes(x = times, y = Effect, group = Edge), 
                color = "#984EA3", linewidth = 0.8, alpha = 0.8, linetype = "dashed") +
      
      # e. 总拟合曲线 (Observed) - 蓝色 #00809D
      geom_line(data = df_total_fit, aes(x = times, y = Total_Fit), 
                color = "#00809D", linewidth = 1.2) +
      
      # g. 末端标签打点
      geom_text_repel(data = df_pw_labels, aes(x = times, y = Effect, label = Edge),
                      color = "#67AE6E", nudge_x = 0.05 * original_range, 
                      direction = "y", hjust = 0, segment.size = 0.2, fontface = "bold", size = 3) +
      
      #geom_text_repel(data = df_ind_labels, aes(x = times, y = Effect, label = Edge),
      #                color = "#DC2525", nudge_x = 0.05 * original_range, 
      #                direction = "y", hjust = 0, segment.size = 0.2, fontface = "bold", size = 3) +
      
      geom_text_repel(data = df_hoi_labels, aes(x = times, y = Effect, label = Edge),
                      color = "#984EA3", nudge_x = 0.05 * original_range, 
                      direction = "y", hjust = 0, segment.size = 0.2, fontface = "bold", size = 3) +
      
      # 【修改 1】：去除 scales = "free_y"，使用默认的 fixed，强制统一坐标轴
      facet_wrap(~ Node, ncol = 4) +
      
      scale_x_continuous(limits = c(min(times), max(times) * 1.001)) + 
      
      theme_bw(base_size = 14) +
      
      theme(
        legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        
        # 【修改 2】：将分面之间的纵向和横向空白间距设置为 0
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "grey80", fill = NA)
        
      ) +
      labs(
        x = "Locus Index",
        y = "Effect"
      )
    
    # print(plot_hoi_2)
    
    
    #ggsave("Active_HOI.pdf", plot_hoi_2,width = 12,height = 8)
    
    return(list(p_high_active = plot_hoi_2,
                df_hoi_active = df_hoi)) 
  }
  
  
  active_HOI = get_active_HOI()
  
  
  ggsave(filename = paste0(name,"decompose.pdf"),active_HOI$p_high_active,width = 12,height = 3)
  
  
  
  f_group2 <- function(jj) {
    
    betam <- beta_all[-c(1:2), jj]
    est <- matrix(0, nrow = length(times_new), ncol = p)
    
    group_idx_list <- split(seq_len(ncol(phi2)), rep(1:p, each = M))
    
    for (g in seq_len(p)) {
      est[, g] <- phi_int2[, group_idx_list[[g]], drop = FALSE] %*% betam[group_idx_list[[g]]]
    }
    
    # 加上截距项和时间线性项
    est[, jj] <- est[, jj] + beta_all[2, jj] * times_new +  beta_all[1, jj]
    
    colnames(est) <- colnames(Y)
    return(est)
  }
  
  f_group_est2 <- lapply(1:p, f_group2)
  
  
  
  f_group <- function(jj) {
    betam <- beta_all[-c(1:2), jj]
    est <- matrix(0, nrow = length(times_norm), ncol = p)
    group_idx_list <- split(seq_len(ncol(phi)), rep(1:p, each = M))
    for (g in seq_len(p)) {
      est[, g] <- phi_int[, group_idx_list[[g]], drop = FALSE] %*% betam[group_idx_list[[g]]]
    }
    est[, jj] <- est[, jj] + beta_all[2, jj] * times_norm +  beta_all[1, jj]
    colnames(est) <- colnames(Y)
    return(est)
  }
  f_group_est <- lapply(1:p, f_group)
  # =========================================================================
  # 1. 提取 Passive 潜在边（同时构建原始网格 edge 与 密级网格 edge2）
  # =========================================================================
  edge <- NULL 
  edge2 <- NULL 
  names = c("H1","H2","S1","S2")
  
  for (s in (1:4)) {
    tmp_mat <- f_group_est[[s]][, -c(s), drop = FALSE] 
    tmp_mat2 <- f_group_est2[[s]][, -c(s), drop = FALSE] 
    
    colnames(tmp_mat) <- paste0(names[s], "<-", colnames(tmp_mat))
    colnames(tmp_mat2) <- paste0(names[s], "<-", colnames(tmp_mat2)) 
    
    if (is.null(edge)) {
      edge <- tmp_mat
      edge2 <- tmp_mat2
    } else {
      edge <- cbind(edge, tmp_mat) 
      edge2 <- cbind(edge2, tmp_mat2) 
    }
  }
  
  col_sums_abs <- colSums(abs(edge))
  edge <- edge[, col_sums_abs != 0, drop = FALSE]
  edge2 <- edge2[, col_sums_abs != 0, drop = FALSE] 
  edge2 <- edge2[,   colMeans(abs(edge2))>effect_thr] 
  
  
  
  # 如果没有任何潜在被动边，直接跳过
  if(ncol(edge) == 0) {
    cat("No passive HOI edges found.\n")
  } else {
    
    # =========================================================================
    # 2. 全局多任务联合矩阵构建 (Simultaneous Estimation Preparation)
    # =========================================================================
    num_passive_edges <- ncol(edge2)
    edge_names <- colnames(edge2)
    
    phi_int_list <- list()    # 原始网格积分基 (用于拟合)
    phi_int2_list <- list()   # 密级网格积分基 (用于绘图)
    y_list <- list()          # 响应变量列表
    group_list <- list()      # ADSIHT 组索引
    candidate_names_list <- list() # 记录每个任务对应的候选变量名
    
    current_group_offset <- 0
    
    for (j in 1:num_passive_edges) {
      main_name = unlist(strsplit(edge_names[j], split = '<-'))[2]
      
      # 2.1 构造当前任务 j 的候选变量 (原始网格)
      candidate_Y = x_smooth[, -match(unlist(strsplit(edge_names[j], split = '<-')), colnames(Y)), drop = FALSE]
      candidate_Y = cbind(x_smooth[, match(main_name, colnames(Y))], candidate_Y)
      colnames(candidate_Y)[1] = main_name
      candidate_names_list[[j]] <- colnames(candidate_Y)
      
      phi_passive <- Reduce(cbind, lapply(1:ncol(candidate_Y), function(xi) {
        poly_basis_1d(candidate_Y[, xi], name = colnames(candidate_Y)[xi], degree = M)
      }))
      phi_passive_int <- apply(phi_passive, 2, function(col) cumtrapz(times_norm, col))
      
      # 2.2 构造当前任务 j 的候选变量 (密级网格)
      candidate_Y2 = x_smooth2[, -match(unlist(strsplit(edge_names[j], split = '<-')), colnames(Y)), drop = FALSE]
      candidate_Y2 = cbind(x_smooth2[, match(main_name, colnames(Y))], candidate_Y2)
      colnames(candidate_Y2)[1] = main_name
      
      phi_passive2 <- Reduce(cbind, lapply(1:ncol(candidate_Y2), function(xi) {
        poly_basis_1d(candidate_Y2[, xi], name = colnames(candidate_Y2)[xi], degree = M)
      }))
      phi_passive_int2 <- apply(phi_passive2, 2, function(col) cumtrapz(times_new, col))
      
      # 2.3 存入全局列表
      phi_int_list[[j]] <- phi_passive_int
      phi_int2_list[[j]] <- phi_passive_int2
      y_list[[j]] <- as.numeric(edge[, j])
      
      # 分配严格递增的 Group ID
      local_groups <- rep(1:ncol(candidate_Y), each = M) + current_group_offset
      group_list[[j]] <- local_groups
      current_group_offset <- max(local_groups)
    }
    
    # 组合为巨大的 Block-Diagonal 设计矩阵
    X_global_bdiag <- as.matrix(bdiag(phi_int_list))
    Y_global_all <- unlist(y_list)
    groups_global <- unlist(group_list)
    
    # =========================================================================
    # 3. 全局联合多任务估计 (Global ADSIHT & Ridge)
    # =========================================================================
    cat(sprintf("Simultaneously estimating %d passive HOI tasks...\n", num_passive_edges))
    fit_global <- ADSIHT(X_global_bdiag, Y_global_all, group = groups_global)
    beta_global <- fit_global$beta[, which.min(fit_global$ic)]
    
    ridge <- TRUE # 是否开启全局 Ridge 微调
    if (ridge) {
      free_idx_global <- which(beta_global != 0)
      if (length(free_idx_global) > 0) {
        X_sub_global <- X_global_bdiag[, free_idx_global, drop = FALSE]
        lambda_global <- 1e-5
        I_global <- diag(length(free_idx_global))
        
        A_global <- t(X_sub_global) %*% X_sub_global + lambda_global * I_global
        b_global <- t(X_sub_global) %*% Y_global_all
        
        beta_ridge_global <- tryCatch({
          solve(A_global, b_global)
        }, error = function(e) {
          MASS::ginv(A_global) %*% b_global
        })
        
        beta_global[free_idx_global] <- as.numeric(beta_ridge_global)
      }
    }
    
    # =========================================================================
    # 4. 参数解包与绘图渲染函数 (按需对各个 target j 进行出图)
    # =========================================================================
    plot_passive_HOI <- function(j) {
      
      # 定位全局 Beta 向量中属于任务 j 的索引范围
      cols_per_task <- sapply(phi_int_list, ncol)
      start_idx <- if(j == 1) 1 else sum(cols_per_task[1:(j-1)]) + 1
      end_idx <- sum(cols_per_task[1:j])
      
      # 提取当前任务的专属参数和基矩阵
      beta_j <- beta_global[start_idx:end_idx]
      phi_int2_j <- phi_int2_list[[j]]
      edge0_obs <- edge[, j]
      c_names <- candidate_names_list[[j]]
      
      # 计算密级网格的总体拟合 (粗红线)
      edge0_hat2 <- phi_int2_j %*% beta_j
      
      # 拆解密级网格上的单节点效应 (彩色虚线)
      node_effects2 <- sapply(1:length(c_names), function(k) {
        idx <- ((k - 1) * M + 1):(k * M)
        return(as.vector(phi_int2_j[, idx, drop = FALSE] %*% beta_j[idx]))
      })
      colnames(node_effects2) <- c_names
      
      # 过滤掉全 0 的无效应节点
      node_effects2 <- node_effects2[, colSums(abs(node_effects2)) != 0, drop = FALSE]
      
      # --- 组装画图数据 (散点分离原则) ---
      t_vec_norm <- as.numeric(times_norm)
      t_vec_new <- as.numeric(times_new)
      
      df_main_obs <- data.frame(times = t_vec_norm, Observed = as.numeric(edge0_obs))
      df_main_fit <- data.frame(times = t_vec_new, Fitted = as.numeric(edge0_hat2))
      
      if(ncol(node_effects2) > 0) {
        df_nodes <- melt(data.frame(times = t_vec_new, node_effects2), 
                         id.vars = "times", variable.name = "Node", value.name = "Effect")
        max_t <- max(df_nodes$times, na.rm = TRUE)
        df_labels <- df_nodes[df_nodes$times == max_t, ]
      } else {
        df_nodes <- data.frame(times = numeric(0), Node = character(0), Effect = numeric(0))
        df_labels <- data.frame()
      }
      
      # --- ggplot 渲染 ---
      p_decomp <- ggplot() +
        geom_line(data = df_nodes, aes(x = times, y = Effect, color = Node), 
                  linewidth = 1, linetype = "dashed", alpha = 0.8) +
        geom_line(data = df_main_fit, aes(x = times, y = Fitted, linetype = "Overall Fitted"), 
                  color = "firebrick", linewidth = 1.5) +
        geom_point(data = df_main_obs, aes(x = times, y = Observed, shape = "Observed Data"), 
                   color = "black", size = 2.5, alpha = 0.6) +
        geom_hline(yintercept = 0, color = "darkgray", linewidth = 0.8) +
        
        geom_text_repel(data = df_labels, 
                        aes(x = times, y = Effect, label = Node, color = Node),
                        nudge_x = 0.05 * diff(range(t_vec_new)), 
                        direction = "y", hjust = 0, segment.size = 0.3,
                        segment.color = "grey50", fontface = "bold", show.legend = FALSE) + 
        
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) + 
        theme_bw(base_size = 14) +
        
        # --- 核心修改区域 ---
        theme(
          legend.position = "none",                   # 1. 彻底隐藏图例
          panel.border = element_blank(),             # 2. 去除四周的方框边框
          axis.line = element_line(color = "black"),  # 3. 补上 L 型的 X 轴和 Y 轴基础实线
          panel.grid.major = element_blank(),         # 4. 去除主要灰色网格线 (选填，使背景全白)
          panel.grid.minor = element_blank()          # 5. 去除次要灰色网格线
        ) +
        # --------------------
      
      labs(
        title = paste("Simultaneous Passive HOI:", edge_names[j]),
        x = "Normalized Time", y = "Effect Size"
        # 移除了 labs 里的 color, shape, linetype，因为图例已经不显示了
      ) +
        scale_shape_manual(values = c("Observed Data" = 16)) +
        scale_linetype_manual(values = c("Overall Fitted" = "solid"))
      
      return(list(p =  p_decomp,
                  data = df_nodes))
    }
  }
  # 批量出图：直接绘制所有任务的图
  passive_cal <- lapply(1:num_passive_edges, function(xi) plot_passive_HOI(xi))
  
  passive_plots = lapply(1:num_passive_edges, function(xi)passive_cal[[xi]][[1]])
  
  
  p_passive_all = wrap_plots(passive_plots,ncol = 4)
  
  passive_HOI_data = lapply(1:num_passive_edges, function(xi)passive_cal[[xi]][[2]])
  names(passive_HOI_data) = colnames(edge2)
  ggsave("passive_HOI.png",p_passive_all,width = 20,height = 7)
  
  
  
  f_group <- function(jj) {
    
    betam <- beta_all[-c(1:2), jj]
    est <- matrix(0, nrow = length(times_new), ncol = p)
    
    group_idx_list <- split(seq_len(ncol(phi2)), rep(1:p, each = M))
    
    for (g in seq_len(p)) {
      est[, g] <- phi_int2[, group_idx_list[[g]], drop = FALSE] %*% betam[group_idx_list[[g]]]
    }
    
    # 加上截距项和时间线性项
    est[, jj] <- est[, jj] + beta_all[2, jj] * times_new +  beta_all[1, jj]
    
    colnames(est) <- colnames(Y)
    return(est)
  }
  
  f_group_est <- lapply(1:p, f_group)
  
  
  adj_main_est <- matrix(0, nrow = p, ncol = p)
  # # 只要该组估计值总和不全为0，即认为存在边(考虑截距，即存在截距认为存在主效应)
  # for (i in 1:p) {
  #   
  #   adj_main_est[i, which(colSums(f_group_est[[i]]) != 0)] <- 1
  # }
  #不考虑截距与初值
  
  # 只要该组beta不全为0，即认为存在边
  for (i in 1:p) {
    tmp = which(rowSums(abs(matrix(beta_all[-c(1,2),i], ncol = M, byrow = TRUE)))!=0)
    adj_main_est[i, tmp] <- 1
  }
  rownames(adj_main_est) = colnames(Y)
  colnames(adj_main_est) = colnames(Y)
  
  
  
  
  
  
  #convert network
  convert_network <- function(xi) {
    effect = as.data.frame(f_group_est[[xi]])
    rownames(effect) = times_new
    ind_name = colnames(effect)[xi]
    
    # ### 新增/修改开始：应用效应阈值过滤 ###
    col_sums = colSums(abs(effect))
    mean_effects = abs(colMeans(effect)) # 计算平均效应的绝对值
    
    # 剔除条件：原来全0的，或者平均效应绝对值小于设定阈值 effect_thr 的
    rm_name = colnames(effect)[col_sums == 0 | mean_effects < effect_thr]
    # ### 新增/修改结束 ###
    
    # 2. 保护机制：无论如何都不能把目标节点自己删掉
    rm_name = setdiff(rm_name, ind_name)
    
    # 3. 安全剔除
    keep_cols = setdiff(colnames(effect), rm_name)
    effect = effect[, keep_cols, drop = FALSE]
    
    # 4. 安全重排
    dep_name = setdiff(keep_cols, ind_name)
    effect = effect[, c(ind_name, dep_name), drop = FALSE]
    
    names(effect)[1] = 'independent'
    effect_mean = colMeans(effect)
    ind_effect = unname(effect_mean[1])
    
    # 5. 安全生成边
    if (length(dep_name) > 0) {
      edge_mean = data.frame(From = dep_name, 
                             To = ind_name, 
                             Effect = effect_mean[-1],
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
    } else {
      edge_mean = data.frame(From = character(0), 
                             To = character(0), 
                             Effect = numeric(0),
                             stringsAsFactors = FALSE,
                             check.names = FALSE)
    }
    return(list(ind_name = ind_name,
                dep_name = dep_name,
                effect_mean = effect_mean,
                edge_mean = edge_mean,
                ind_effect = ind_effect))
  }
  
  network_data = lapply(1:p,convert_network)
  
  active_HOI_data = active_HOI$df_hoi_active
  active_HOI_data$Edge = gsub("..", "<-", active_HOI_data$Edge, fixed = TRUE)
  
  
  active_HOI_sum <- aggregate(Effect ~ Edge + Node, 
                              data = active_HOI_data, 
                              FUN = function(x) mean(x, na.rm = TRUE))
  
  colnames(active_HOI_sum) = c("From","To","Effect")
  
  
  
  
  
  summary_passive <- function(d,name){
    if (length(unique(d$Node)) == 1) {
      return(NULL)
    } else{
      current_name = substr(name, 5, 6)
      df_summary_base <- aggregate(Effect ~ Node, data = d, 
                                   FUN = function(x) mean(x, na.rm = TRUE))
      df_summary_base = df_summary_base[-match(current_name,df_summary_base$Node),]
      # 2. 规范列名
      names(df_summary_base)[1] <- "From"
      df_summary_base$To =   name
      return(df_summary_base)
    }
  }
  
  
  passive_HOI_sum = lapply(1:ncol(edge2),function(xi) 
    summary_passive( passive_HOI_data[[xi]],name = colnames(edge2)[xi]))
  passive_HOI_sum = Reduce(rbind,passive_HOI_sum)
  
  
  
  network_data = list(pairwise = network_data,
                      active = active_HOI_sum,
                      pasive = passive_HOI_sum)
  
  return(list(beta_all = beta_all,
              times = times,
              times_new = times_restored,
              Y= Y,
              f_group_est = f_group_est, 
              network_data = network_data,
              adj_est = adj_main_est,
              fig0 = pp0,
              fig = pp1))
}





Y = t(sapply(1:18,function(xi) rowMeans(module.df[[xi]])))



colnames(Y) = c("H1","H2","S1","S2")

times = rowSums(Y)
new_times = seq(min(times),max(times),length = 18)
times = new_times
smooth = "bs"
M=3
effect_thr = 1e-3
name = "Overall"

res = MTODE(Y,new_times,M =3, smooth = "bs",name = "Overall")
result = res$network_data

network_plot_HOI <- function(result, title = NULL, maxeffect = NULL){
  require(igraph)
  
  #plot mean network
  
  normalization <- function(x, z = 0.1){(x-min(x))/(max(x)-min(x))+z}
  
  extra <- sapply(result$pairwise,"[[", "ind_effect")
  after <- data.frame(do.call(rbind, lapply(result$pairwise, "[[", "edge_mean")))
  
  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL
  
  
  after$edge.colour <- ifelse(after$Effect >= 0, "#D50000", "#0033CC")
  
  active = result$active[,c(1,3,2)]
  active$edge.colour <- ifelse(active$Effect >= 0, "#D50000", "#0033CC")
  
  
  passive = result$pasive[,c(1,3,2)]
  passive$edge.colour <- ifelse(passive$Effect >= 0, "#D50000", "#0033CC")
  
  
  after$Type = "pairwise"
  active$Type = "active"
  passive$Type = "passive"
  
  
  all_edges = rbind( after,active,passive )
  
  
  
  #nodes
  nodes <- data.frame(sapply(result$pairwise, '[[','ind_name'),
                      sapply(result$pairwise, '[[','ind_name'),
                      extra)
  
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence = 0.1
  nodes$influence[match(aggregate(Effect ~ To, data = after, sum)[,1],nodes$id)] <- 
    aggregate(Effect ~ To, data = after, sum)[,2]
  
  nodes$node.colour = NA
  for (i in 1:nrow(nodes)) {
    if(nodes$influence[i]>=0){
      nodes$node.colour[i] = "#F9F6C4"
    } else{
      nodes$node.colour[i] = "#F9F6C4"
    }
  }
  
  
  
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  
  #layout
  n_nodes <- vcount(net)
  l <- matrix(0, nrow = n_nodes, ncol = 2)
  
  # 修改 1：定义长方形的基础顶点 (拉长 X 轴，Y 轴保持不变，形成 2:1 的长方形)
  base_pos <- data.frame(
    name = c("S1", "S2", "H1", "H2"),
    x = c(-1.2, 1.2, -1.2, 1.2),
    y = c(1, -1, -1, 1)
  )
  
  theta <- -20 * pi / 180  
  
  for (i in 1:n_nodes) {
    node_name <- V(net)$name[i]
    
    idx <- match(node_name, base_pos$name)
    if (!is.na(idx)) {
      x <- base_pos$x[idx]
      y <- base_pos$y[idx]
    } else {
      x <- 0; y <- 0 
    }
    
    # 剪切变换
    x_shear <- x + 0.5 * y
    y_shear <- y
    
    # 旋转变换
    x_rot <- x_shear * cos(theta) - y_shear * sin(theta)
    y_rot <- x_shear * sin(theta) + y_shear * cos(theta)
    
    l[i, 1] <- x_rot
    l[i, 2] <- y_rot
  }
  
  # 修改 2：动态计算绘图边界，加上 0.5 的余量防止节点被切边
  x_margin <- c(min(l[,1]) - 0.5, max(l[,1]) + 0.5)
  y_margin <- c(min(l[,2]) - 0.5, max(l[,2]) + 0.5)
  
  
  
  
  
  
  base_nodes <- nodes$name
  hyper_nodes <- union(all_edges$From,all_edges$To)[grep("<",union(all_edges$From,all_edges$To))]
  
  
  # 创建空图并添加所有顶点
  g <- make_empty_graph(directed = TRUE) + 
    vertices(c(base_nodes, hyper_nodes))
  
  
  
  
  # ==========================================
  # 1. 重构图对象与基础属性
  # ==========================================
  all_vertices <- data.frame(name = c(base_nodes, hyper_nodes), stringsAsFactors = FALSE)
  g <- graph_from_data_frame(d = all_edges, vertices = all_vertices, directed = TRUE)
  
  # 【色彩恢复】：直接使用你原始计算好的红蓝色彩映射
  E(g)$color <- E(g)$edge.colour
  E(g)$Type <- all_edges$Type
  
  # 【区分主动与被动】：通过线型 (lty) 来清晰区分高阶调控的类型
  # 1 = Solid(实线, 基础边), 2 = Dashed(虚线, 主动调控), 4 = Dot-dash(点划线, 被动反馈)
  E(g)$lty <- ifelse(E(g)$Type == "pairwise", 1, 
                     ifelse(E(g)$Type == "active", 2, 4))
  
  # 设定弯曲度：只让 Pairwise 基础边弯曲，高阶调控的连线保持笔直
  # (注意：正数表示逆时针弯曲，负数表示顺时针弯曲)
  base_curvature <- 0.1 
  E(g)$curved <- ifelse(E(g)$Type == "pairwise", base_curvature, 0.05)
  
  # 粗细映射
  abs_effect <- abs(as.numeric(E(g)$Effect))
  E(g)$width <- (abs_effect / max(abs_effect + 1e-6)) * 6 + 2
  E(g)$arrow.size <- 1.0
  
  # ==========================================
  # 2. 核心数学：计算 3D 贝塞尔曲线的绝对中点
  # ==========================================
  layout_full <- matrix(NA, nrow = vcount(g), ncol = 2)
  rownames(layout_full) <- V(g)$name
  
  # 2.1 填入基础节点 3D 坐标 (来自你的透视底板 l)
  for (i in 1:vcount(net)) {
    node_name <- V(net)$name[i]
    layout_full[node_name, ] <- l[i, ]
  }
  
  # 2.2 计算超边节点在曲线上的精确位置
  for (h_node in hyper_nodes) {
    # 拆解出源节点和目标节点 (例如 "S1<-H2" -> 目标 S1, 源 H2)
    parts <- unlist(strsplit(h_node, "<|-")) 
    target <- parts[1]
    source <- parts[3] 
    
    p1 <- layout_full[source, ]
    p2 <- layout_full[target, ]
    
    # 获取方向向量
    dx <- p2[1] - p1[1]
    dy <- p2[2] - p1[2]
    
    # igraph 二次贝塞尔曲线的中点偏移数学公式
    offset_x <- -dy * base_curvature * 0.2
    offset_y <-  dx * base_curvature * 0.2
    
    # 完美钉在弧线的最顶端！
    layout_full[h_node, 1] <- (p1[1] + p2[1])/2 + offset_x
    layout_full[h_node, 2] <- (p1[2] + p2[2])/2 + offset_y
  }
  
  # ==========================================
  # 3. 映射顶点属性 (隐藏超边节点)
  # ==========================================
  is_hyper <- grepl("<", V(g)$name)
  
  # 【核心修改】：将超边节点的形状设为 "none" (无)，大小设为 0
  V(g)$shape <- ifelse(is_hyper, "none", "circle")
  V(g)$size <- ifelse(is_hyper, 0, 30) 
  
  # 【核心修改】：将超边节点的颜色和边框彻底设为 NA (透明)
  base_color_map <- setNames(nodes$node.colour, nodes$name)
  V(g)$color <- ifelse(is_hyper, NA, base_color_map[V(g)$name])
  V(g)$frame.color <- ifelse(is_hyper, NA, "black")
  
  # 标签处理：超边依然不显示标签
  V(g)$label <- ifelse(is_hyper, "", V(g)$name)
  V(g)$label.color <- "black"
  V(g)$label.font <- 2
  V(g)$label.cex <- ifelse(grepl("<", V(g)$name), 0, 1.4)
  
  # ==========================================
  # 4. 终极渲染：三子网左右排列
  # ==========================================
  # 将颜色变浅的辅助函数（混入白色）
  lighten_color <- function(col, factor = 0.72) {
    rgb_val <- col2rgb(col) / 255
    rgb_light <- rgb_val + (1 - rgb_val) * factor
    rgb(rgb_light[1,], rgb_light[2,], rgb_light[3,])
  }
  
  plot_subnet <- function(graph, sub_lty, main_title) {
    # 选出目标类型的边
    keep_edges <- which(E(graph)$lty == sub_lty)
    supplemented_pw <- integer(0)  # 记录补入的 pairwise 边 id
    
    if (sub_lty != 1) {
      # 找出这些边涉及的超节点（名称含 "<"）
      ep <- ends(graph, keep_edges, names = TRUE)
      involved_nodes <- union(ep[, 1], ep[, 2])
      hyper_involved <- involved_nodes[grepl("<", involved_nodes)]
      
      # 对每个超节点补入对应的 pairwise 边
      if (length(hyper_involved) > 0) {
        pw_edges <- which(E(graph)$lty == 1)
        pw_ep    <- ends(graph, pw_edges, names = TRUE)
        for (hn in hyper_involved) {
          parts  <- unlist(strsplit(hn, "<|-"))
          target <- parts[1]
          source <- parts[3]
          matched <- pw_edges[pw_ep[, 1] == source & pw_ep[, 2] == target]
          supplemented_pw <- union(supplemented_pw, matched)
        }
        keep_edges <- union(keep_edges, supplemented_pw)
      }
    }
    
    sg <- subgraph.edges(graph, eids = keep_edges, delete.vertices = FALSE)
    sg_names <- V(sg)$name
    lf <- layout_full[sg_names, , drop = FALSE]
    
    edge_colors <- E(sg)$color
    edge_lty    <- E(sg)$lty
    
    if (sub_lty != 1) {
      for (i in seq_len(ecount(sg))) {
        ep_i <- ends(sg, i, names = TRUE)
        is_hyper_edge <- grepl("<", ep_i[1]) || grepl("<", ep_i[2])
        
        if (E(sg)$lty[i] == 1) {
          # 补入的 pairwise 基础边：弱化颜色
          edge_colors[i] <- lighten_color(E(sg)$color[i])
        } else if (is_hyper_edge) {
          # 超边（HOI 调控连线）：保持原色，改为实线
          edge_lty[i] <- 1
        }
      }
    }
    
    xm <- c(min(lf[,1]) - 0.5, max(lf[,1]) + 0.5)
    ym <- c(min(lf[,2]) - 0.5, max(lf[,2]) + 0.5)
    
    plot(sg,
         layout             = lf,
         xlim               = xm,
         ylim               = ym,
         rescale            = FALSE,
         main               = main_title,
         vertex.size        = V(sg)$size * 3 + 0.5,
         vertex.label.color = "black",
         edge.curved        = E(sg)$curved,
         edge.color         = edge_colors,
         edge.lty           = edge_lty)
  }
  
  pdf(paste0(title, '.pdf'), width = 12, height = 4)
  par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
  plot_subnet(g, 1, "Pairwise (Solid)")
  plot_subnet(g, 2, "Active HOI (Dashed)")
  plot_subnet(g, 4, "Passive HOI (Dot-dash)")
  par(mfrow = c(1, 1))
  dev.off()
  
  
  
  ###################################
  
  node_names <- V(g)$name
  node_types <- ifelse(grepl("<", node_names), "HyperEdge", "BaseNode")
  
  
  df_nodes <- data.frame(
    Node = node_names,
    Type = node_types,
    
    # 度中心性 (Degree)：连接数
    Degree_Total = degree(g, mode = "all"),
    Degree_In    = degree(g, mode = "in"),  # 被调控的数量
    Degree_Out   = degree(g, mode = "out"), # 调控别人的数量
    
    # 介数中心性 (Betweenness)：充当最短路径“桥梁”的次数
    Betweenness = betweenness(g, directed = TRUE, weights = NA, normalized = TRUE),
    
    # 接近中心性 (Closeness)：到其他所有节点的平均距离的倒数
    Closeness = closeness(g, mode = "all", weights = NA, normalized = TRUE),
    
    # PageRank (网页排名算法)：评估节点的重要程度（被重要节点指向则更重要）
    PageRank = page_rank(g, directed = TRUE, weights = NA)$vector
  )
  
  # 将结果按 PageRank 或 Betweenness 降序排列，找出网络中的“核心枢纽”
  df_nodes <- df_nodes[order(df_nodes$PageRank, decreasing = TRUE), ]
  rownames(df_nodes) <- NULL
  
  cat("\n--- 节点级别指标 (前 10 名核心节点) ---\n")
  print(head(df_nodes, 10))
  
  
  df_graph <- data.frame(
    Metric = c(
      "网络密度 (Edge Density)", 
      "网络直径 (Diameter)", 
      "平均最短路径 (Average Path Length)", 
      "全局聚类系数 (Global Transitivity)"
    ),
    Value = c(
      edge_density(g),                                   # 实际边数与可能最大边数的比例
      diameter(g, directed = TRUE, weights = NA),        # 网络中最长的最短路径
      mean_distance(g, directed = TRUE, weights = NA),   # 所有节点对之间的平均最短路径
      transitivity(g, type = "global")                   # 形成三角形(反馈环)的概率
    )
  )
  
  cat("\n--- 全局网络指标 ---\n")
  print(df_graph)
  
  
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  cat("\n开始绘制网络拓扑学特征分析图...\n")
  
  # 统一的颜色主题 (BaseNode 用亮黄/橙，HyperEdge 用深灰/蓝)
  color_palette <- c("BaseNode" = "#FF9800", "HyperEdge" = "#607D8B")
  
  # ==========================================
  # 图 A: Top 15 核心节点排名 (棒棒糖图)
  # ==========================================
  # 提取前 15 名节点
  top_nodes <- df_nodes %>% head(15)
  
  p1 <- ggplot(top_nodes, aes(x = reorder(Node, PageRank), y = PageRank, color = Type)) +
    geom_segment(aes(xend = reorder(Node, PageRank), y = 0, yend = PageRank), 
                 color = "grey80", size = 1) +
    geom_point(size = 4) +
    coord_flip() +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    labs(title = "A. Top Nodes by PageRank", 
         x = "Nodes (Base & Hyper)", 
         y = "PageRank Score") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none",
          panel.grid.major.y = element_blank())
  
  # ==========================================
  # 图 B: 基础节点 vs 超边节点 的度分布 (箱线图)
  # ==========================================
  # 整理数据为长格式，方便画分面图
  df_deg <- df_nodes %>%
    select(Node, Type, Degree_In, Degree_Out) %>%
    pivot_longer(cols = c(Degree_In, Degree_Out), 
                 names_to = "Degree_Type", 
                 values_to = "Value")
  
  # 清洗标签名
  df_deg$Degree_Type <- ifelse(df_deg$Degree_Type == "Degree_In", "In-Degree", "Out-Degree")
  
  p2 <- ggplot(df_deg, aes(x = Type, y = Value, fill = Type, color = Type)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.5) +
    facet_wrap(~Degree_Type, scales = "free_y") +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    theme_bw() +
    labs(title = "B. Degree Distribution",
         x = "", y = "Number of Connections") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none",
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"))
  
  # ==========================================
  # 图 C: 网络中心性景观图 (散点图)
  # ==========================================
  # 挑选需要打标签的极值节点 (PageRank极高 或 Betweenness极高)
  highlight_nodes <- df_nodes %>%
    filter(PageRank > quantile(PageRank, 0.9) | Betweenness > quantile(Betweenness, 0.9))
  
  p3 <- ggplot(df_nodes, aes(x = Betweenness, y = PageRank, color = Type, size = Degree_Total)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(data = highlight_nodes, aes(label = Node), 
                    size = 3.5, fontface = "bold", color = "black",
                    box.padding = 0.5, max.overlaps = 15, show.legend = FALSE) +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(range = c(2, 8)) +
    theme_classic() +
    labs(title = "C. Centrality Landscape",
         x = "Betweenness Centrality", 
         y = "PageRank",
         size = "Total Degree") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  
  # ==========================================
  # 4. 使用 Patchwork 拼图并输出
  # ==========================================
  # 布局设计：上面两张图并排，下面一张大图横跨
  final_plot <- (p1 | p2) / p3 + 
    plot_annotation(
      title = "Topology Analysis of High-Order Interaction Network",
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  
  
  return(final_plot)
}



network_plot_HOI(res$network_data, title = "Overall")
