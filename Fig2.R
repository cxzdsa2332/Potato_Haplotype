rm(list = ls())

library(MASS)
library(reshape2)
library(ADSIHT)
library(ggplot2)
library(MASS)
library(orthopolynom)
library(pracma)
library(Matrix)
library(ggrepel)
library(patchwork)

df = read.csv(file = "log2_FPKM_3_replicates_at_4_haplotypes.csv", row.names = 1)
df = df[,-c(1,2)]
df = df-min(df)


haplotype_groups <- list(grep("H1",colnames(df)), 
                         grep("H2",colnames(df)), 
                         grep("S1",colnames(df)), 
                         grep("S2",colnames(df)))


dfs = list(H1 = df[,haplotype_groups[[1]]],
           H2 = df[,haplotype_groups[[2]]],
           S1 = df[,haplotype_groups[[3]]],
           S2 = df[,haplotype_groups[[4]]])
dfs2 = data.frame(t(Reduce(cbind,lapply(dfs, rowSums))))
rownames(dfs2) = names(dfs)
dfs2 <- dfs2[, colSums(dfs2, na.rm = TRUE) != 0]



power_equation <- function(x, dat_par){ t(sapply(1:nrow(dat_par),function(c) dat_par[c,1]*x^dat_par[c,2] ) )}


power_equation_base <- function(y, times) {
  x <- as.numeric(times)
  y <- as.numeric(y)
  
  # ---------------------------------------------------------
  # 第一步：获取高质量的初始值 (去除 0 的影响)
  # ---------------------------------------------------------
  # 仅提取 y > 0 且 x > 0 的有效点用于初始线性拟合
  valid_idx <- (y > 0) & (x > 0)
  
  # 如果有效点太少，直接返回 NULL 避免报错
  if (sum(valid_idx) < 2) {
    warning("有效非零数据点少于2个，无法拟合")
    return(NULL)
  }
  
  x_valid <- x[valid_idx]
  y_valid <- y[valid_idx]
  
  # 用干净的数据进行 log-log 线性拟合
  lmFit <- lm(log(y_valid) ~ log(x_valid))
  coefs <- coef(lmFit)
  
  # 提取初始值
  a_start <- exp(coefs[1])
  b_start <- coefs[2]
  
  # ---------------------------------------------------------
  # 第二步：使用非线性最小二乘法进行最终拟合
  # ---------------------------------------------------------
  # 将数据整理为 data.frame，这是 nls 的最佳实践
  df <- data.frame(x = x, y = y)
  
  # 运行带边界约束的 nls
  model <- try(
    nls(
      y ~ a * (x^b), 
      data = df,
      start = list(a = as.numeric(a_start), b = as.numeric(b_start)),
      lower = c(a = 1e-10, b = -Inf), # 约束 a 必须为正数
      algorithm = "port",             # 启用边界约束的算法
      control = nls.control(
        maxiter = 1000, 
        minFactor = 1e-10, 
        warnOnly = TRUE               # 遇到轻微不收敛时返回最优结果而不是直接报错
      )
    ), 
    silent = TRUE
  )
  
  # 返回结果
  if (inherits(model, "try-error")) {
    return(NULL)
  } else {
    return(model)
  }
}

times = colSums(dfs2)

fit = lapply(1:4,function(xi) power_equation_base(dfs2[xi,],times))
pars = t(sapply(fit, function(xi) coef(xi)))
df_hat = power_equation(times,pars)
rownames(df_hat) = rownames(dfs2)


# ── 绘图 ──────────────────────────────────────────────────────────────────────

# 颜色方案：点统一色，拟合线区分 S/H
col_point <- c(S1 = "#58CAF1", S2 = "#58CAF1", H1 = "#58CAF1", H2 = "#58CAF1")
col_line  <- c(S1 = "#03829E", S2 = "#03829E", H1 = "#03829E", H2 = "#03829E")

# 构建散点数据（observed）
scatter_df <- do.call(rbind, lapply(rownames(dfs2), function(hap) {
  data.frame(
    haplotype = hap,
    x         = as.numeric(times),
    y         = as.numeric(dfs2[hap, ]),
    stringsAsFactors = FALSE
  )
}))

# 构建拟合曲线数据（smooth power-law line）
x_seq <- seq(min(times), max(times), length.out = 400)
line_df <- do.call(rbind, lapply(seq_len(nrow(pars)), function(i) {
  hap <- rownames(dfs2)[i]
  a   <- pars[i, 1]
  b   <- pars[i, 2]
  data.frame(
    haplotype = hap,
    x         = x_seq,
    y_fit     = a * x_seq^b,
    stringsAsFactors = FALSE
  )
}))

# 固定 facet 顺序
hap_order <- c("S1", "S2", "H1", "H2")
scatter_df$haplotype <- factor(scatter_df$haplotype, levels = hap_order)
line_df$haplotype    <- factor(line_df$haplotype,    levels = hap_order)

p <- ggplot() +
  geom_point(
    data    = scatter_df,
    aes(x = x, y = y, color = haplotype),
    size    = 1.5,
    alpha   = 0.35,
    shape   = 16
  ) +
  geom_line(
    data      = line_df,
    aes(x = x, y = y_fit, color = haplotype),
    linewidth = 1.3
  ) +
  scale_color_manual(
    values = c(col_point, col_line)[hap_order],
    breaks = hap_order
  ) +
  ggnewscale::new_scale_color() +
  geom_line(
    data      = line_df,
    aes(x = x, y = y_fit, color = haplotype),
    linewidth = 1.3
  ) +
  facet_wrap(~ haplotype, ncol = 4) +
  labs(
    x = "Locus Index",
    y = "Effect"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "none",
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12),
    panel.spacing     = unit(0.3, "lines"),
    panel.border      = element_rect(color     = "grey60",
                                     fill      = NA,
                                     linewidth = 0.5),
    axis.line         = element_blank()
  )

# 点和线分别用独立颜色：重新组织为两层独立 aes
.p1_bg <- data.frame(
  haplotype = factor(hap_order, levels = hap_order),
  bg        = c("#FEF5F6", "#FEF5F6", "#F2F9FC", "#F2F9FC"),
  stringsAsFactors = FALSE
)

p1 <- ggplot() +
  geom_rect(
    data = .p1_bg,
    aes(fill = bg),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_point(
    data  = scatter_df,
    aes(x = x, y = y, color = haplotype),
    size  = 1.5,
    alpha = 0.6,
    shape = 16
  ) +
  scale_color_manual(values = col_point) +
  ggnewscale::new_scale_color() +
  geom_line(
    data      = line_df,
    aes(x = x, y = y_fit, color = haplotype),
    linewidth = 1.3
  ) +
  scale_color_manual(values = col_line) +
  facet_wrap(~ haplotype, ncol = 4) +
  labs(
    x = "Locus Index",
    y = "Effect"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "none",
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12),
    panel.spacing.x   = unit(0, "lines"),
    panel.spacing.y   = unit(0.3, "lines"),
    panel.border      = element_rect(color = "grey60", fill = NA, linewidth = 0.5),
    axis.line         = element_blank()
  )






legendre_basis_1d <- function(x, name = NULL, degree = 3) {
  # x 必须在 [0, 1] 之间
  n <- length(x)
  P <- matrix(0, n, degree)
  
  # 移位 Legendre 多项式公式 (Shifted Legendre Polynomials on [0,1])
  # P_1(x) = 2x - 1
  # P_2(x) = 6x^2 - 6x + 1
  # P_3(x) = 20x^3 - 30x^2 + 12x - 1
  # 递推公式: (n+1)P_{n+1}(x) = (2n+1)(2x-1)P_n(x) - nP_{n-1}(x)
  
  # P0 = 1 (通常包含在截距中，这里从 1阶开始)
  p_curr <- 2 * x - 1       # P1
  p_prev <- rep(1, n)       # P0
  
  P[, 1] <- p_curr
  
  if (degree >= 2) {
    for (d in 2:degree) {
      k <- d - 1
      # 递推公式计算 P_{k+1} 即当前的 d
      p_next <- ((2 * k + 1) * (2 * x - 1) * p_curr - k * p_prev) / (k + 1)
      P[, d] <- p_next
      
      p_prev <- p_curr
      p_curr <- p_next
    }
  }
  
  base_name <- if (!is.null(name)) name else "x"
  colnames(P) <- paste0(base_name, "_Leg_", 1:degree)
  return(P)
}



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





#p = 4
Y <- as.matrix(t(dfs2))
times = rowSums(Y)
Y = Y[order(times),]
times = rowSums(Y)
M = 3
smooth = "power_equation"



MTODE <- function(Y, times, M = 5, smooth = "bs",effect_thr = 1e-1) {
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
  lambda <- 1e-5
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
        
        # If you want to distinguish which 's' this comes from, you can modify the column names, for example:
        colnames(tmp_mat) <- paste0(names[s], "<-", colnames(tmp_mat))
        
        if (is.null(edge)) {
          edge <- tmp_mat
        } else {
          edge <- cbind(edge, tmp_mat) # Concatenate horizontally by columns, automatically retaining column names
        }
      }
      col_sums_abs <- colSums(abs(edge))
      edge <- edge[, col_sums_abs != 0, drop = FALSE]
      return(edge)
    }
    
    edges <- lapply(1:4, function(xi) get_possible_edge(xi))
    
    
    edges <- lapply(edges, function(df) {
      col_means <- colMeans(abs(df))
      
      keep_cols <- names(col_means)[col_means >= effect_thr*10]
      
      df_filtered <- df[, keep_cols, drop = FALSE]
      
      return(df_filtered)
    })
    
    
    num_blocks <- length(edges) 
    
    # 1. Iterate through the 4 elements, storing the generated phi_high_int in a list
    phi_high_int_list <- lapply(1:num_blocks, function(i) {
      
      # Generate polynomial basis for each column of the i-th edge matrix
      phi_list_high <- lapply(1:ncol(edges[[i]]), function(xi) {
        poly_basis_1d(edges[[i]][, xi], name = colnames(edges[[i]])[xi], degree = M)
      })
      
      # Merge by column
      phi_high <- Reduce(cbind, phi_list_high)
      
      # Integrate the basis
      phi_high_int <- apply(phi_high, 2, function(col) cumtrapz(times_norm, col))
      
      return(phi_high_int)
    })
    
    # 2. Assemble the list of 4 matrices into a block diagonal matrix
    # bdiag() returns a sparse matrix (dgCMatrix) by default; use as.matrix() to convert it to a standard dense matrix for subsequent routine operations
    phi_high_int_bdiag <- as.matrix(bdiag(phi_high_int_list))
    high_group <- rep(1:ncol(Reduce(cbind, edges)), each = M)
    
    # Extract global error and perform joint ADSIHT selection
    error <- as.vector(residuals_mat)
    fit_high <- ADSIHT(phi_high_int_bdiag, error, group = high_group)
    beta_high <- fit_high$beta[, which.min(fit_high$ic)]
    
    
    
    free_idx_high <- which(beta_high != 0)
    beta_high_final <- beta_high
    
    if (length(free_idx_high) > 0) {
      
      # Extract the valid feature subset (drop = FALSE is mandatory to prevent demotion error for a single column)
      X_sub_high <- phi_high_int_bdiag[, free_idx_high, drop = FALSE]
      y_sub_high <- error
      
      lambda_high <- 1e-150
      
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
        colnames(dep_mat) <- paste0(colnames(Y)[-i], "->", node_name)
        df_pw <- melt(data.frame(times = times, dep_mat), id.vars = "times",
                      variable.name = "Edge", value.name = "Effect")
        df_pw$Node <- node_name
        df_pairwise_list[[i]] <- df_pw
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
    
    # ── 过滤小效应────────────────────────────────────────────────
    .p3_thr <- 0.1
    # pairwise：按 Edge×Node 组计算均值绝对值，移除低于阈值的组
    if (nrow(df_pairwise) > 0) {
      .pw_mean <- tapply(abs(df_pairwise$Effect),
                         interaction(df_pairwise$Node, df_pairwise$Edge),
                         mean, na.rm = TRUE)
      .pw_keep <- names(.pw_mean)[.pw_mean >= .p3_thr]
      .pw_key  <- interaction(df_pairwise$Node, df_pairwise$Edge)
      df_pairwise <- df_pairwise[.pw_key %in% .pw_keep, ]
    }
    # HOI：同理
    if (nrow(df_hoi) > 0) {
      .hoi_mean <- tapply(abs(df_hoi$Effect),
                          interaction(df_hoi$Node, df_hoi$Edge),
                          mean, na.rm = TRUE)
      .hoi_keep <- names(.hoi_mean)[.hoi_mean >= .p3_thr]
      .hoi_key  <- interaction(df_hoi$Node, df_hoi$Edge)
      df_hoi    <- df_hoi[.hoi_key %in% .hoi_keep, ]
    }
    
    # ==========================================
    # 提取末端标签数据
    # ==========================================
    max_t_orig <- max(times, na.rm = TRUE)
    
    df_pw_labels <- df_pairwise[df_pairwise$times == max_t_orig, ]
    
    if (nrow(df_hoi) > 0) {
      max_t_hoi <- max(df_hoi$times, na.rm = TRUE)
      df_hoi_labels <- df_hoi[abs(df_hoi$times - max_t_hoi) < 1e-6, ]
    } else {
      df_hoi_labels <- data.frame(times = numeric(0), Edge = character(0), Effect = numeric(0), Node = character(0))
    }
    
    # 截断/替换标签文字
    df_pw_labels$Edge <- substr(as.character(df_pw_labels$Edge), 1, 2)
    if (nrow(df_hoi_labels) > 0)
      df_hoi_labels$Edge <- gsub("\\.\\.", "<", as.character(df_hoi_labels$Edge))
    
    # ==========================================
    # 绘制最终的 plot_hoi_2
    # ==========================================
    # 统一 y 轴范围
    .p3_all_y <- c(
      df_obs$Observed, df_pairwise$Effect,
      df_ind$Effect,   df_total_fit$Total_Fit
    )
    if (nrow(df_hoi) > 0) .p3_all_y <- c(.p3_all_y, df_hoi$Effect)
    .p3_ymin <- floor(min(.p3_all_y, na.rm = TRUE))
    .p3_ymax <- ceiling(max(.p3_all_y, na.rm = TRUE))
    
    # 分面背景色
    .p3_bg <- data.frame(
      Node = factor(target_order, levels = target_order),
      bg   = c("#FEF3F3", "#FEF3F3", "#F0F7FE", "#F0F7FE"),
      stringsAsFactors = FALSE
    )
    
    plot_hoi_2 <- ggplot() +
      # 0. 分面背景色
      geom_rect(
        data = .p3_bg,
        aes(fill = bg),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
        inherit.aes = FALSE
      ) +
      scale_fill_identity() +
      # f. 原始观测数据
      geom_point(data = df_obs, aes(x = times, y = Observed),
                 color = "#53CBF3", size = 1.5, alpha = 0.4) +
      # a. 零线基准
      geom_hline(yintercept = 0, linetype = 2, color = "grey60",
                 linewidth = 0.4) +
      # b. Pairwise 效应 - 绿色
      geom_line(data = df_pairwise,
                aes(x = times, y = Effect, group = Edge),
                color = "#67AE6E", linewidth = 0.8, alpha = 0.8) +
      # c. Independent 效应 - 红色
      geom_line(data = df_ind,
                aes(x = times, y = Effect, group = Edge),
                color = "#DC2525", linewidth = 1) +
      # d. HOI 效应 - 紫色虚线
      geom_line(data = df_hoi,
                aes(x = times, y = Effect, group = Edge),
                color = "#984EA3", linewidth = 0.55, alpha = 0.9,
                linetype = "longdash") +
      # e. 总拟合曲线 - 蓝色
      geom_line(data = df_total_fit,
                aes(x = times, y = Total_Fit),
                color = "#00809D", linewidth = 1.2) +
      # g. 末端标签
      geom_text_repel(
        data = df_pw_labels, aes(x = times, y = Effect, label = Edge),
        color = "#67AE6E", nudge_x = 0.05 * original_range,
        direction = "y", hjust = 0, segment.size = 0.2,
        fontface = "bold", size = 3) +
      geom_text_repel(
        data = df_hoi_labels, aes(x = times, y = Effect, label = Edge),
        color = "#984EA3", nudge_x = 0.05 * original_range,
        direction = "y", hjust = 0, segment.size = 0.2,
        fontface = "bold", size = 3) +
      facet_wrap(~ Node, ncol = 4) +
      scale_x_continuous(
        limits = c(min(times) * 0.96, max(times) * 1.15)) +
      scale_y_continuous(limits = c(.p3_ymin, .p3_ymax)) +
      labs(x = "Locus Index", y = "Effect") +
      theme_bw(base_size = 13) +
      theme(
        legend.position  = "none",
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        panel.border     = element_rect(
          color = "grey55", fill = NA,
          linewidth = 0.45),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 12,
                                        color = "black"),
        panel.spacing.x  = unit(0, "lines"),
        panel.spacing.y  = unit(0.3, "lines"),
        axis.line        = element_blank()
      )
    
    #ggsave("Active_HOI.pdf", plot_hoi_2, width = 16, height = 5)
    return(plot_hoi_2)
    
  }
  
  
  
  active_HOI = get_active_HOI()
  
  
  
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
      
      # --- 组装画图数据：使用原始 Locus Index ---
      t_vec_obs <- as.numeric(times)
      t_vec_new <- as.numeric(times_restored)
      
      df_main_obs <- data.frame(times = t_vec_obs, Observed = as.numeric(edge0_obs))
      
      if(ncol(node_effects2) > 0) {
        df_nodes <- melt(data.frame(times = t_vec_new, node_effects2),
                         id.vars = "times", variable.name = "Node", value.name = "Effect")
        # 按整体均值决定颜色：正→紫红，负→蓝紫
        .nd_mean <- tapply(df_nodes$Effect, as.character(df_nodes$Node), mean, na.rm = TRUE)
        .nd_col  <- ifelse(.nd_mean >= 0, "#B5338A", "#5B5EA6")
        df_nodes$line_col <- .nd_col[as.character(df_nodes$Node)]
        
        max_t <- max(df_nodes$times, na.rm = TRUE)
        df_labels <- df_nodes[df_nodes$times == max_t, ]
      } else {
        df_nodes <- data.frame(times = numeric(0), Node = character(0),
                               Effect = numeric(0), line_col = character(0))
        df_labels <- data.frame()
      }
      
      # --- ggplot 渲染 ---
      p_decomp <- ggplot() +
        geom_line(data = df_nodes, aes(x = times, y = Effect,
                                       group = Node, color = line_col),
                  linewidth = 1, alpha = 0.85) +
        geom_point(data = df_main_obs, aes(x = times, y = Observed),
                   color = "#6aaed6", size = 2, alpha = 0.5) +
        geom_hline(yintercept = 0, color = "darkgray", linewidth = 0.5,
                   linetype = 2) +
        
        geom_text_repel(data = df_labels,
                        aes(x = times, y = Effect, label = Node,
                            color = line_col),
                        nudge_x = 0.05 * diff(range(t_vec_new)),
                        direction = "y", hjust = 0, segment.size = 0.3,
                        segment.color = "grey50", fontface = "bold",
                        show.legend = FALSE) +
        scale_color_identity() +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
        theme_bw(base_size = 14) +
        theme(
          legend.position  = "none",
          panel.border     = element_blank(),
          axis.line        = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        labs(
          title = paste("Passive HOI:", edge_names[j]),
          x = "Locus Index", y = "Effect Size"
        )
      
      
    }
  }
  # 批量出图：只选取第2和第5个任务
  passive_plots <- lapply(1:num_passive_edges, plot_passive_HOI)
  .sel_idx <- intersect(c(3, 5), seq_along(passive_plots))
  #.sel_idx = 1:7
  p_passive_all = wrap_plots(passive_plots[.sel_idx], ncol = 2)
  
  #ggsave("passive_HOI.png",p_passive_all,width = 20,height = 13)
  
  
  
  
  
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
  
  
  return(list(beta_all = beta_all,
              times = times,
              times_new = times_restored,
              Y= Y,
              f_group_est = f_group_est,
              network_data = network_data,
              adj_est = adj_main_est,
              fig0 = pp0,
              fig = pp1,
              plot_hoi_2 = active_HOI,
              p_passive = p_passive_all))
}

plot_decompose <- function(res, j, values = c('#67AE6E','#DC2525','#00809D'), effect_thr = 1e-1) {
  require(ggrepel)
  require(reshape2) # 确保 melt 可用
  require(ggplot2)
  
  orig_times = res$times
  plot_df = data.frame(res$f_group_est[[j]])
  
  # 修复潜在的外部环境变量引用问题：使用 res$Y 而不是 Y
  ind_name = colnames(res$Y)[j] 
  
  # ### 过滤太低的效应 ###
  # 在计算 obs(总和) 之前，先计算各列真实效应
  mean_effects = abs(colMeans(plot_df))
  col_sums = colSums(abs(plot_df))
  
  # 找出全0或平均效应绝对值低于阈值的变量
  rm_name = colnames(plot_df)[col_sums == 0 | mean_effects < effect_thr]
  # 保护机制：永远不要把目标节点自身删掉
  rm_name = setdiff(rm_name, ind_name)
  
  rownames(plot_df) = res$times_new
  
  # 计算总估计值 obs (此时还包含所有变量，保证拟合总线的准确性)
  plot_df$obs = rowSums(plot_df) 
  
  # 实施剔除
  if (length(rm_name) > 0) {
    keep_cols = setdiff(colnames(plot_df), rm_name)
    plot_df = plot_df[, keep_cols, drop = FALSE]
  }
  
  # 转换数据格式
  plot_df2 = melt(as.matrix(plot_df))
  plot_df2$Var2 = as.character(plot_df2$Var2) # 强转字符，防止因子类型导致报错
  
  plot_df2$type = 'dependent'
  plot_df2$type[plot_df2$Var2 == ind_name] = 'independent'
  plot_df2$type[plot_df2$Var2 == 'obs'] = 'observed'
  
  # 构建标签数据框 (获取最后一行的时间点数据用于打标签)
  label_df = data.frame(Var1 = as.numeric(plot_df[nrow(plot_df), ]))
  rownames(label_df) = colnames(plot_df)
  label_df$name = rownames(label_df)
  
  # 移除 obs 和自身 (ind_name)，只给其他的 dependent 打标签
  label_df = label_df[label_df$name != 'obs', , drop = FALSE]
  label_df = label_df[label_df$name != ind_name, , drop = FALSE]
  
  # 设置标签的 X 轴位置 (排在最右侧)
  label_df$Var2 = max(orig_times)
  
  # 开始绘图
  pp = ggplot() + 
    geom_point(mapping = aes(x = orig_times, y = res$Y[, j]), color = 'blue', alpha = 0.3) +
    geom_line(plot_df2, mapping = aes(x = Var1, y = value, group = Var2, color = type), linetype = 1) +
    theme_bw() + 
    ylab('Values') + 
    ggtitle(paste0(ind_name, "_dynamics")) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(limits = c(min(orig_times) * 0.96, max(orig_times) * 1.15)) +
    scale_color_manual(values = values) +
    
    # ### 新增：去除图例和背景网格 ###
    theme(
      legend.position = "none",           # 去除图例
      panel.grid.major = element_blank(), # 去除主网格线
      panel.grid.minor = element_blank()  # 去除次网格线
    )
  
  # 只有存在需要打标签的变量时，才添加 geom_text_repel，防止 label_df 为空时报错
  if (nrow(label_df) > 0) {
    pp = pp + geom_text_repel(data = label_df,
                              aes(x = Var2,
                                  y = Var1,
                                  label = name),
                              color = values[1], # 对应 dependent 的颜色
                              size = 3.5,
                              nudge_x = max(orig_times) * 0.1,
                              max.overlaps = Inf)
  }
  
  return(pp)
}
res = MTODE(Y,times,M = 3, smooth = "power_equation")



# ── p2：统一分面版，舍弃 patchwork ────────────────────────────────────────────
hap_indices <- c(3, 4, 1, 2)          # S1, S2, H1, H2
hap_order   <- c("S1", "S2", "H1", "H2")
.effect_thr <- 1e-1

lines_list  <- list()
points_list <- list()
labels_list <- list()

for (.idx in seq_along(hap_indices)) {
  .j        <- hap_indices[.idx]
  .hap      <- hap_order[.idx]
  .ind_name <- colnames(res$Y)[.j]
  .orig_t   <- res$times
  
  .pdf <- data.frame(res$f_group_est[[.j]])
  .meff <- abs(colMeans(.pdf))
  .csums <- colSums(abs(.pdf))
  .rm <- setdiff(colnames(.pdf)[.csums == 0 | .meff < .effect_thr], .ind_name)
  
  rownames(.pdf) <- res$times_new
  .pdf$obs <- rowSums(.pdf)
  if (length(.rm) > 0)
    .pdf <- .pdf[, setdiff(colnames(.pdf), .rm), drop = FALSE]
  
  .m <- reshape2::melt(as.matrix(.pdf))
  .m$Var1 <- as.numeric(as.character(.m$Var1))
  .m$Var2 <- as.character(.m$Var2)
  .m$type <- ifelse(.m$Var2 == .ind_name, "independent",
                    ifelse(.m$Var2 == "obs", "observed", "dependent"))
  .m$haplotype <- .hap
  lines_list[[.idx]] <- .m
  
  points_list[[.idx]] <- data.frame(
    x = .orig_t, y = as.numeric(res$Y[, .j]), haplotype = .hap)
  
  .ldf <- data.frame(y_end = as.numeric(.pdf[nrow(.pdf), ]),
                     name  = colnames(.pdf), stringsAsFactors = FALSE)
  .ldf <- .ldf[!.ldf$name %in% c("obs", .ind_name), , drop = FALSE]
  if (nrow(.ldf) > 0) {
    .ldf$x <- max(.orig_t); .ldf$haplotype <- .hap
    labels_list[[.idx]] <- .ldf
  }
}

.all_lines  <- do.call(rbind, lines_list)
.all_points <- do.call(rbind, points_list)
.all_labels <- if (length(labels_list) > 0) do.call(rbind, labels_list) else NULL

# 统一因子顺序
.all_lines$haplotype  <- factor(.all_lines$haplotype,  levels = hap_order)
.all_points$haplotype <- factor(.all_points$haplotype, levels = hap_order)
if (!is.null(.all_labels))
  .all_labels$haplotype <- factor(.all_labels$haplotype, levels = hap_order)

# dependent 线按整体均值着色：正→紫红，负→蓝紫
.dep_rows <- .all_lines$type == "dependent"
.dep_mean <- tapply(.all_lines$value[.dep_rows],
                    .all_lines$Var2[.dep_rows],
                    mean, na.rm = TRUE)
.dep_col  <- ifelse(.dep_mean >= 0, "#B5338A", "#5B5EA6")
.all_lines$line_col <- ifelse(
  .dep_rows,
  .dep_col[as.character(.all_lines$Var2)],
  ifelse(.all_lines$type == "independent", "#DC2525", "#00809D")
)

# 统一 y 轴：最低取数据最小值，最高固定 90
.y_min <- floor(min(c(.all_lines$value, .all_points$y), na.rm = TRUE))
.y_max <- 90

# 每个分面的背景色
.bg_df <- data.frame(
  haplotype = factor(hap_order, levels = hap_order),
  bg        = c("#FEF3F3", "#FEF3F3", "#F0F7FE", "#F0F7FE"),
  stringsAsFactors = FALSE
)

p2 <- ggplot() +
  # 1. 分面背景色矩形
  geom_rect(
    data = .bg_df,
    aes(fill = bg),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  # 2. 散点（原始观测）
  geom_point(
    data  = .all_points,
    aes(x = x, y = y),
    color = "#6aaed6", alpha = 0.35, size = 1.2
  ) +
  # 3. 零线
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.4) +
  # 4. 效应线
  geom_line(
    data = .all_lines,
    aes(x = Var1, y = value, group = Var2, color = line_col),
    linewidth = 0.9
  ) +
  scale_color_identity() +
  # 5. 末端标签
  {if (!is.null(.all_labels) && nrow(.all_labels) > 0)
    geom_text_repel(
      data  = .all_labels,
      aes(x = x, y = y_end, label = name),
      color = "#4a9a58", size = 3.2,
      nudge_x    = diff(range(res$times)) * 0.08,
      direction  = "y", hjust = 0,
      segment.size = 0.25, max.overlaps = Inf
    )
  } +
  # 6. 分面与轴
  facet_wrap(~ haplotype, ncol = 4) +
  scale_x_continuous(
    limits = c(min(res$times) * 0.96, max(res$times) * 1.14)
  ) +
  scale_y_continuous(limits = c(.y_min, .y_max)) +
  labs(x = NULL, y = "Values") +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "none",
    panel.background  = element_blank(),
    panel.grid        = element_blank(),
    panel.border      = element_rect(
      color = "grey55", fill = NA, linewidth = 0.45),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12),
    panel.spacing.x   = unit(0, "lines"),
    panel.spacing.y   = unit(0.3, "lines"),
    axis.line         = element_blank()
  )


p3 <- res$plot_hoi_2

p4 = res$p_passive




network_plot <- function(result, title = NULL, maxeffect = NULL, dd = 1.2){
  require(igraph)
  #plot mean network
  
  normalization <- function(x, z = 0.1){(x-min(x))/(max(x)-min(x))+z}
  
  extra <- sapply(result,"[[", "ind_effect")
  after <- data.frame(do.call(rbind, lapply(result, "[[", "edge_mean")))
  
  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL
  
  
  after$edge.colour <- ifelse(after$Effect >= 0, "#D50000", "#0033CC")
  #nodes
  nodes <- data.frame(sapply(result, '[[','ind_name'),
                      sapply(result, '[[','ind_name'),
                      extra)
  
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence[match(aggregate(Effect ~ To, data = after, sum)[,1],nodes$id)] <- 
    aggregate(Effect ~ To, data = after, sum)[,2]
  
  # 所有节点→优美淡黄色
  nodes$node.colour <- "#FAF0BE"
  
  #normalization
  if (is.null(maxeffect)) {
    after[,3] <- normalization(abs(after[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  } else{
    after[,3] <- (abs(after[,3]))/maxeffect*1+0.5
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1+0.5
  }
  
  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )
  
  #layout
  n_nodes <- vcount(net)
  l <- matrix(0, nrow = n_nodes, ncol = 2)
  
  # 节点位置：左上S1，右上H2，右下S2，左下H1；dd控制左右间距
  base_pos <- data.frame(
    name = c("S1", "H2", "S2", "H1"),
    x    = c(-dd,  dd,  dd, -dd),
    y    = c(  1,   1,  -1,  -1)
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
  
  plot.igraph(net,
              rescale = FALSE,               # 修改 3：关闭自动缩放，保留长方形比例
              xlim = x_margin,               # 传入计算好的 X 边界
              ylim = y_margin,               # 传入计算好的 Y 边界
              vertex.label = V(net)$name,
              vertex.label.color = "black",
              vertex.shape = "circle",
              vertex.label.cex = V(net)$ind_effect*1.5 + 0.5,
              vertex.size = V(net)$ind_effect * 30 + 50,
              edge.arrow.size = E(net)$Effect*2 + 0.8,
              edge.curved = 0.05,
              edge.color = E(net)$edge.colour,
              edge.frame.color = E(net)$edge.colour,
              edge.width = E(net)$Effect*5+1,
              vertex.color = V(net)$node.colour,
              layout = l, 
              main = title
  )
  return(net)
}

pdf(file = "Fig2_network_plot.pdf", width = 8, height = 6)

# 2. 执行绘图代码 (此时图表不会显示在屏幕上，而是直接画入 PDF)
network_plot(res$network_data)

# 3. 关闭图形设备，完成文件保存 (非常重要，漏掉这句文件会损坏)
dev.off()
ggsave('Fig2A.pdf',p1,width = 13,height=3)
