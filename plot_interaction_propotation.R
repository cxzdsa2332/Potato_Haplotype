rm(list = ls())
load("funclu.RData")

library(igraph)
library(pracma)
library(ADSIHT)
library(Matrix)
library(reshape2)
library(ggplot2)




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







MTODE0 <- function(Y, times, M = 3, smooth = "bs") {
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
    
  } else if(smooth == "power_equation"){
    
    # 使用 sapply 遍历每一列进行非线性拟合
    x_smooth <- sapply(1:p, function(xi) {
      y_val <- as.numeric(Y[, xi])
      
      # 加上微小偏移量防止 0 的负指数导致无穷大
      x_val <- times_norm + 1e-6 
      
      # 使用 tryCatch 处理拟合不收敛的情况
      tryCatch({
        # 1. 初始值估计 (Initial Guesses)
        # 假设 y = a * x^b
        # a 估计为范围，b 初始设为 1 (线性)
        c_start <- min(y_val)
        a_start <- max(y_val) - min(y_val)
        if(a_start == 0) a_start <- 0.1 # 防止常数数列导致错误
        
        # 2. 非线性最小二乘拟合 (NLS)
        fit_nls <- nls(y_val ~ a * I(x_val^b), 
                       start = list(a = a_start, b = 1),
                       control = nls.control(maxiter = 100, warnOnly = TRUE))
        
        # 3. 返回预测值
        predict(fit_nls)
        
      }, error = function(e) {
        # 4. 回退机制 (Fallback)
        # 如果数据不符合幂律分布（例如是波动的），NLS 会失败。
        # 此时回退到局部加权回归 (LOESS) 以保证程序不中断。
        warning(paste("Var", colnames(Y)[xi], ": Power fit failed, using Loess instead."))
        fitted(loess(y_val ~ times_norm, span = 0.75))
      })
    })
    # 保持列名一致
    colnames(x_smooth) <- colnames(Y)
  }
  
  
  # 2. 构建基函数并积分
  # 构造多项式基
  phi_list <- lapply(1:p, function(xi) {
    poly_basis_1d(x_smooth[, xi], name = colnames(Y)[xi], degree = M)
  })
  phi <- Reduce(cbind, phi_list)
  
  # 对基进行积分 (Trapezoidal integration)
  phi_int <- apply(phi, 2, function(col) cumtrapz(times_norm, col))
  
  
  
  
  
  cal_Locus <- function(i){
    main = colnames(Y)[i]
    dep =  colnames(Y)[-i]
    p <- ncol(Y)
    y = Y[,i]
    # --- 1. 准备设计矩阵与分组 ---
    # 确保 phi_int 是矩阵
    X_basis <- as.matrix(phi_int) 
    # 合并时间项作为协变量（通常作为控制变量）
    X_design <- cbind(Time = as.numeric(times_norm), X_basis)
    
    # 更新分组向量：假设 Time 项自成一组（组号 1），其余基函数按原逻辑分组
    # 注意：很多 R 包要求 group 从 1 开始且为连续整数
    terms_per_var <- sapply(phi_list, ncol)
    groups_vec <- c(1, rep(2:(p + 1), times = terms_per_var))
    
    # --- 2. 数据标准化 (Standardization) ---
    # 使用 scale 虽方便，但手动计算均值和标准差在还原时更稳健
    X_center <- colMeans(X_design)
    X_scale  <- apply(X_design, 2, sd)
    X_scale[X_scale == 0] <- 1 # 防止除以 0
    
    X_design_scaled <- scale(X_design, center = X_center, scale = X_scale)
    
    y_vec <- as.numeric(y)
    y_center <- mean(y_vec)
    y_scale  <- sd(y_vec)
    if(y_scale == 0) y_scale <- 1
    y_scaled <- (y_vec - y_center) / y_scale
    
    # --- 3. 调用 ADSIHT ---
    fit <- ADSIHT(X_design_scaled, y_scaled, group = groups_vec,s0 = 1)
    
    # 选择最佳模型 (基于 IC 信息准则)
    # 注意：需检查 fit$ic 是否存在，部分版本可能叫 fit$BIC 或 fit$AIC
    if(!is.null(fit$ic)){
      best_idx <- which.min(fit$ic)
    } else {
      #  # 备选方案：取最后一个 lambda 或指定的稀疏度
      best_idx <- ncol(fit$beta) 
    }
    
    # 提取标准化后的系数
    coef_scaled <- fit$beta[, best_idx]
    # --- 4. 还原系数 (De-normalization) ---
    # 1. 还原斜率： beta_orig = beta_scaled * (sd_y / sd_x)
    coef_restored <- coef_scaled * (y_scale / X_scale)
    
    # 2. 还原截距： beta0 = mean_y - sum(mean_x * beta_orig)
    intrinsic_intercept_scaled <- if(!is.null(fit$intercept)) fit$intercept[best_idx] else 0
    intercept_restored <- (intrinsic_intercept_scaled * y_scale + y_center) - sum(X_center * coef_restored)
    
    #alternative by sgl
    #library(sparsegl)
    #fit_sgl_cv <- cv.sparsegl(x = X_design_scaled, 
    #                          y = y_scaled, 
    #                          group = groups_vec,
    #                          asparse = 1)
    #coef_scaled <- as.matrix(coef(fit_sgl_cv, s = "lambda.1se"))
    #
    #coef_restored <- coef_scaled[-1] * (y_scale / X_scale)
    
    # 2. 还原截距： beta0 = mean_y - sum(mean_x * beta_orig)
    # 特别注意：如果 ADSIHT 内部有自己的截距项(fit$intercept)，需要加上去
    #intrinsic_intercept_scaled <- coef_scaled[1]
    #intercept_restored <- (intrinsic_intercept_scaled * y_scale + y_center) - sum(X_center * coef_restored)
    
    
    
    # --- 5. 结果整理 ---
    final_coefs <- c(Intercept = intercept_restored, coef_restored)
    
    
    
    idx_exclude <- grep(main, colnames(X_design))
    X_sub <- X_design[, -c(1,idx_exclude)]
    coef_sub <- coef_restored[-c(1,idx_exclude)]
    components_matrix <- sweep(X_sub, 2, coef_sub, "*")
    
    
    col_labels <- rep(colnames(Y)[-i], times = M)
    compressed_matrix <- sapply(colnames(Y)[-i], function(var_name) {
      sub_cols <- components_matrix[, grep(var_name,col_labels)]
      rowSums(sub_cols)
    })
    
    
    idx <- which(colSums(compressed_matrix) != 0)
    
    if ( length(idx) ==0 ) {
      pp = NULL
      find_type = NULL
    } else{
      
      # 加上 drop = FALSE 防止仅有一个变量时退化为向量并丢失列名
      dep_value_mat <- compressed_matrix[, idx, drop = FALSE] 
      
      # --- 计算拟合值与独立分量 ---
      y_fitted <- intercept_restored + as.numeric(X_design %*% coef_restored)
      
      # 计算 ind_value
      idx_main <- grep(main, colnames(X_design))
      ind_value <- intercept_restored + 
        X_design[, 1] * coef_restored[1] + 
        as.numeric(X_design[, idx_main, drop = FALSE] %*% coef_restored[idx_main])
      
      
      
      find_dep = colnames(dep_value_mat)
      
      main_locus = strsplit(main,'_')[[1]][2]
      dep_locus  = sapply(strsplit(find_dep,'_'),'[[', 2)
      
      
      
      
      
      find_type <- rep(NA,length(dep_locus))
      if ( any(main_locus == dep_locus) ){
        #如果有相同位点，先考虑相同位点
        tmp = which(dep_locus == main_locus)
        current_main = substr(main,1,1)
        current_dep =  substr(find_dep[tmp],1,1)
        find_type0 <- ifelse(current_dep == current_main, 
                             "Auto Dominance", 
                             "Allo Dominance")
        
        
        
        find_type[tmp] = find_type0
        
        dep_locus = dep_locus[-tmp]
        
        
        if (length(dep_locus)>0) {
          current_dep = substr(find_dep[-tmp],1,1)
          find_type0 <- ifelse(current_dep == current_main, 
                               "Cis Epistasis", 
                               "Auto Trans Epistasis")
          find_type[-tmp] = find_type0
          
          
          
        } else{
          
        }
        
      } else{
        
        current_main = substr(main,1,1)
        current_dep =  substr(find_dep,1,1)
        
        
        find_type <- ifelse(current_dep == current_main, 
                            "Auto Dominance", 
                            "Allo Dominance")
        
      }
      
      
      idx <- which(find_type == "Cis Epistasis")
      
      # 如果存在至少一个 "Cis Epistasis"
      if (length(idx) > 0) {
        
        main_prefix <- substr(main, 1, 2)
        
        # 提取对应的 find_dep 的首字母 (这会是一个向量)
        dep_prefixes <- substr(find_dep[idx], 1, 2)
        
        # 使用 ifelse 进行向量化条件判断和批量赋值
        find_type[idx] <- ifelse(main_prefix == dep_prefixes,
                                 "Cis Epistasis",
                                 "Allo Trans Epistasis")
      }
      
      
      # --- 准备绘图数据 (Base Data) ---
      library(ggplot2)
      
      plot_data <- data.frame(
        Time = times, 
        Original = as.numeric(y),
        Fitted = y_fitted,
        IndValue = ind_value
      )
      
      # --- 构建多 dep_value 的长数据格式 (用于多线绘制) ---
      dep_df <- data.frame(
        Time = rep(times, times = ncol(dep_value_mat)),
        Variable = rep(colnames(dep_value_mat), each = length(times)),
        Value = as.vector(dep_value_mat)
      )
      
      # --- 准备标签数据 (取每条 dep_value 线的最后一个点) ---
      last_time <- max(times)
      label_data <- data.frame(
        Time = rep(last_time, ncol(dep_value_mat)),
        Value = dep_value_mat[nrow(dep_value_mat), ],
        Label = colnames(dep_value_mat)
      )
      
      # --- 绘制定制化对比图 ---
      # 计算 x 轴的扩展范围，防止末端标签被截断 (扩展约 10%)
      x_range <- max(times) - min(times)
      x_limit_max <- max(times) + x_range * 0.15 
      
      pp <- ggplot() +
        # 0. 增加 y=0 的辅助虚线 (建议放在底层，以免遮挡数据线)
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.8) +
        # 1. 原始数据 (蓝色散点)
        geom_point(data = plot_data, aes(x = Time, y = Original, color = "Original"), 
                   size = 2, alpha = 0.5) +
        # 2. 拟合曲线 (蓝色粗实线)
        geom_line(data = plot_data, aes(x = Time, y = Fitted, color = "Fitted"), 
                  linewidth = 1.2) +
        # 3. 独立变量贡献 ind_value (红色实线)
        geom_line(data = plot_data, aes(x = Time, y = IndValue, color = "Ind_Value"), 
                  linewidth = 1.2) +
        # 4. 依赖变量贡献 dep_value (多条绿色细线以示区分)
        geom_line(data = dep_df, aes(x = Time, y = Value, group = Variable, color = "Dep_Value"), 
                  linewidth = 0.8, alpha = 0.8) +
        # 5. 在 dep_value 线条末端加上文本标签
        geom_text(data = label_data, aes(x = Time, y = Value, label = Label), 
                  color = "#6FAF4F", hjust = -0.2, vjust = 0.5, size = 4, fontface = "bold") +
        # 统一管理颜色映射 (虽然不显示图例，但这部分仍用来控制线条颜色)
        scale_color_manual(
          name = "Legend",
          values = c(
            "Original" = "blue",      # 原始数据蓝
            "Fitted" = "blue",        # 拟合线蓝
            "Ind_Value" = "red",      # ind_value 红
            "Dep_Value" ="#6FAF4F"    # dep_value 绿
          )
        ) +
        # 扩展 X 轴给标签留出空间
        coord_cartesian(xlim = c(min(times), x_limit_max), clip = "off") +
        labs(
          x = "Locus Index",
          y = "Effect"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none", 
          plot.margin = margin(t = 10, r = 20, b = 10, l = 10),
          
          # <--- 新增：去除灰色背景格子(网格线) --->
          panel.grid.major = element_blank(), # 去除主网格线
          panel.grid.minor = element_blank(), # 去除次网格线
          panel.background = element_rect(fill = "white", color = NA) # 确保背景纯白
        )
      
      
    }  
    
    
    
    
    obj <- list(p=pp,
                find_type = find_type,
                main = main,
                n = length(which(colSums(compressed_matrix) != 0)))
    
    obj 
  }
  
  
  # cal_Locus(137)
  
  res = lapply(1:ncol(Y), function(xi) cal_Locus(xi))
  
  
  
  
  interaction_summary = table(unlist(sapply(res,'[[',2)))
  
  
  
  
  
  
  plot_net1 <- function(){
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  return(interaction_summary)
}



get_interacgtion <- function(select){
  #select = 6
  
  
  
  df = result$cluster[result$cluster$apply.omega..1..which.max. == select,]
  df = df[,-ncol(df)]
  
  
  H1 = df[,1:3]
  H2 = df[,4:6]
  S1 = df[,7:9]
  S2 = df[,10:12]
  
  construct_L <- function(xi){
    L1 = rbind(c(H1[xi,]),c(H2[xi,]),c(S1[xi,]),c(S2[xi,]))
    colnames(L1) = c("A","B","C")
    rownames(L1) = c("H1","H2","S1","S2")
    return(L1)
  }
  
  dfs = lapply(1:nrow(H1), construct_L)
  names(dfs) = rownames(df)
  
  
  
  power_equation <- function(x, dat_par){ t(sapply(1:nrow(dat_par),function(c) dat_par[c,1]*x^dat_par[c,2] ) )}
  
  
  power_equation_base <- function(y, times, offset = 1e-6) {
    x <- as.numeric(times)
    y <- as.numeric(y)
    
    # ---------------------------------------------------------
    # 预处理：防止 x <= 0 导致的算术运算崩溃 (Inf/NaN)
    # ---------------------------------------------------------
    # 如果 x 中存在 <= 0 的值，平移整个 x 轴，保证底数绝对大于0
    if (any(x <= 0, na.rm = TRUE)) {
      x <- x - min(x, na.rm = TRUE) + offset
    }
    
    # ---------------------------------------------------------
    # 第一步：获取高质量的初始值 (去除 0 和 NA 的影响)
    # ---------------------------------------------------------
    valid_idx <- (y > 0) & (x > 0) & !is.na(x) & !is.na(y)
    
    # 幂函数至少需要2个参数，建议至少需要3个点以保留自由度
    if (sum(valid_idx) < 3) { 
      return(NULL) 
    }
    
    x_valid <- x[valid_idx]
    y_valid <- y[valid_idx]
    
    # 用干净的数据进行 log-log 线性拟合
    # 使用 suppressWarnings 防止完全共线性时的警告刷屏
    lmFit <- suppressWarnings(lm(log(y_valid) ~ log(x_valid)))
    coefs <- coef(lmFit)
    
    # 容错机制：如果 lm 失败产生 NA，提供备用初值
    if (any(is.na(coefs)) || any(is.infinite(coefs))) {
      a_start <- mean(y_valid, na.rm = TRUE)
      b_start <- 1.0 # 默认退化为线性
    } else {
      a_start <- exp(coefs[1])
      b_start <- coefs[2]
    }
    
    # 防止初值溢出导致 nls 启动失败
    if (a_start == 0) a_start <- 1e-5
    if (is.infinite(a_start)) a_start <- max(y_valid)
    
    # ---------------------------------------------------------
    # 第二步：使用非线性最小二乘法进行最终拟合
    # ---------------------------------------------------------
    df <- data.frame(x = x, y = y)
    
    # 使用 tryCatch 进行更安全的错误拦截
    model <- tryCatch({
      nls(
        y ~ a * (x^b), 
        data = df,
        start = list(a = as.numeric(a_start), b = as.numeric(b_start)),
        # 适当放大 a 的下界。1e-30 在某些机器精度下等同于 0，可能导致奇异梯度
        lower = c(a = 1e-50, b = -Inf), 
        upper = c(a = 100, b = Inf), 
        #algorithm = "port",              
        control = nls.control(
          maxiter = 1000, 
          minFactor = 1e-30, 
          warnOnly = TRUE              
        )
      )
    }, error = function(e) {
      # 捕获到错误直接返回 NULL，不中断后续程序
      return(NULL)
    })
    
    return(model)
  }
  
  
  
  
  times = log2(colSums(2^H1)+colSums(2^H2)+colSums(2^S1)+colSums(2^S2)+1)
  
  
  fit_batch <- function(cc){
    ys = cbind(as.numeric(dfs[[cc]][1,]),
               as.numeric(dfs[[cc]][2,]),
               as.numeric(dfs[[cc]][3,]),
               as.numeric(dfs[[cc]][4,]))
    fit = lapply(1:4,function(xi) power_equation_base(ys[,xi],times))
    
    if (any(sapply(fit, is.null)) == TRUE){
      fitted = matrix(0,nrow = 30,ncol = 4)
      
      no = which(sapply(fit, is.null)!=TRUE)
      fitted0 = lapply(no,function(xi)   
        predict(fit[[xi]], newdata = data.frame(x = seq(min(times, na.rm = TRUE), 
                                                        max(times, na.rm = TRUE), 
                                                        length.out = 30))))
      fitted0 = Reduce(cbind, fitted0)
      fitted[,no] = fitted0
      rownames(fitted) = NULL
      
    } else{
      fitted = lapply(1:4,function(xi)   
        predict(fit[[xi]], newdata = data.frame(x = seq(min(times, na.rm = TRUE), 
                                                        max(times, na.rm = TRUE), 
                                                        length.out = 30))))
      fitted = Reduce(cbind, fitted)
      rownames(fitted) = NULL
    }
    return(fitted)
  }
  
  
  
  fitted_bathch = lapply(1:length(dfs), function(xi) fit_batch(xi))
  for (i in 1:length(dfs)) {
    colnames(fitted_bathch[[i]]) = paste0(c("H1","H2","S1","S2"),"_",names(dfs)[i])
  }
  
  Y = Reduce(cbind,fitted_bathch) #should be Module Locus number*4
  Y = Y[,-which(colSums(Y)==0)] #remove 0
  new_times = seq(min(times),max(times), length = 30)
  
  
  
  interaction = MTODE0(Y,new_times)
  possible_type = c("Allo Dominance","Allo Trans Epistasis","Auto Dominance",
                    "Cis Epistasis","Auto Trans Epistasis")
  
  interaction_df = data.frame(possible_type,num = 0)
  idx = match(names(interaction),possible_type)
  interaction_df$num[idx] = interaction
  
  
  interaction_df
}

#Y = Y;times = new_times;smooth = "bs";M = 3

res = lapply(c(6,7),function(xi) get_interacgtion(xi))



plot_data <- do.call(cbind, res)
plot_data  = plot_data[,-seq(3,ncol(plot_data),by = 2)]
colnames(plot_data)[2:ncol(plot_data)] = paste0("M",1:(ncol(plot_data)-1))




models <- colnames(plot_data)[-1]

# 利用 R 矩阵按列拉平的特性，构建长数据框
plot_data_long <- data.frame(
  Type = rep(plot_data[, 1], times = length(models)),
  Model = rep(models, each = nrow(plot_data)),
  Count = as.vector(as.matrix(plot_data[, -1]))
)

library(ggplot2)

p <- ggplot(plot_data_long, aes(x = Model, y = Count, fill = Type)) +
  
  # position = "fill" 会自动把绝对数值转换为 0-1 的比例进行堆叠
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  
  # 纯 Base R 函数自定义标签，无需依赖 scales 包即可显示 %
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  
  scale_fill_brewer(palette = "Set2") +
  
  labs(
    title = "Interaction Type Proportion per Model",
    x = "Model",
    y = "Proportion",
    fill = "Interaction Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "right"
  )

# 显示图形
print(p)


