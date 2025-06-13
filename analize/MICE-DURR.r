# 必要パッケージ
library(glmnet)

# カテゴリ変数の最頻値関数
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# MICE-DURR アルゴリズム実装
mice_durr <- function(data, M = 5, maxit = 10, alpha = 1) {
  # data: 欠測値を含む data.frame（数値／因子混在可）
  # M: 作成する多重代入データセット数
  # maxit: MICE 内部の反復回数
  # alpha: glmnet の elastic net 混合パラメータ (1=LASSO, 0=Ridge)
  
  n <- nrow(data)
  p <- ncol(data)
  imp_list <- vector("list", M)
  
  for (m in seq_len(M)) {
    # ——— 初期値で初期化 ———
    imp <- data
    for (j in seq_len(p)) {
      miss_j <- is.na(imp[[j]])
      if (is.numeric(imp[[j]])) {
        imp[miss_j, j] <- mean(imp[[j]], na.rm = TRUE)
      } else {
        imp[miss_j, j] <- Mode(imp[[j]])
      }
    }
    
    # ——— MICE イテレーション ———
    for (iter in seq_len(maxit)) {
      for (j in seq_len(p)) {
        # 変数 j の欠測／完全データマスク
        obs_idx <- which(!is.na(data[[j]]))
        mis_idx <- which(is.na(data[[j]]))
        if (length(mis_idx) == 0) next
        
        # W_j: 説明変数行列（z_j を除く全変数、最新の imputed 値を含む）
        W <- as.matrix(imp[, -j, drop = FALSE])
        
        # (1) ブートストラップサンプルを作成
        boot_idx <- sample(obs_idx, length(obs_idx), replace = TRUE)
        X_boot  <- W[boot_idx, , drop = FALSE]
        y_boot  <- imp[boot_idx, j]
        
        # (2) glmnet でモデルフィッティング → θ^(m)
        cvfit <- cv.glmnet(X_boot, y_boot, alpha = alpha)
        # 最適 λ
        lam <- cvfit$lambda.min
        
        # 残差分散を推定
        y_hat_obs <- predict(cvfit, W[obs_idx, , drop = FALSE], s = lam)
        sigma <- sqrt(mean((imp[obs_idx, j] - y_hat_obs)^2))
        
        # (3) 欠測箇所に予測＋ノイズを追加
        X_mis <- W[mis_idx, , drop = FALSE]
        y_pred <- as.numeric(predict(cvfit, X_mis, s = lam))
        imp[mis_idx, j] <- y_pred + rnorm(length(mis_idx), 0, sigma)
      }
    }
    
    imp_list[[m]] <- imp
  }
  
  names(imp_list) <- paste0("Imputation_", seq_len(M))
  return(imp_list)
}

# ——— 使い方例 ———
# ——  MICE-DURR を実行 ——  
# 欠測のある z1–z3 を含む全変数を imputate  
imp_list <- mice_durr(sim_data,
                      M     = 5,    # 多重データセット数
                      maxit = 10,   # 各チェーン内の反復回数
                      alpha = 1)    # 1=LASSO, 0=Ridge

# —— 3. 結果確認例 ——  
# 1つ目の補完データセットの先頭行
head(imp_list[[1]])