#*n   <- 100 #サンプル数
#*p   <- 200 #変数総数  
#*rho <- 0.1 #自己相関
#*q   <- 4   #活性変数数

# —— 0. セットアップ ——  
library(MASS)      # mvrnorm 用  
set.seed(20250611)  # 再現性確保  

# —— 1. パラメータ ——  
n   <- 100         # サンプル数  
p   <- 200         # 変数総数  
rho <- 0.1         # AR(1) 自己相関 

# —— 2. z4–zp の生成 ——  
#   Σ_{ij} = rho^{|i-j|} の (p-3)x(p-3) 共分散行列を作成  
Sigma  <- rho^abs(outer(1:(p-3), 1:(p-3), "-"))  
Z_rest <- mvrnorm(n, mu = rep(0, p-3), Sigma = Sigma)  
colnames(Z_rest) <- paste0("z", 4:p)  

# —— 3. z1–z3 の生成 ——  
S4        <- c(4, 5, 50, 51)     # 真の活性変数番号  
alpha4    <- rep(1, length(S4))  # α = (1,1,1,1)'  
# Z_rest の列番号は z_k → (k-3)  
Z_S       <- Z_rest[, S4 - 3]    # n×4 行列  
# 平均 μ = 1 + Z_S %*% α, 分散 = 4 → sd = 2  
mu_z123   <- 1 + Z_S %*% alpha4  
z123      <- matrix(
  rnorm(n * 3, mean = as.vector(mu_z123), sd = 2),
  nrow = n, ncol = 3
)  
colnames(z123) <- paste0("z", 1:3)  

# —— 4. y の生成 ——  
# y = β0 + ∑_{i=1}^5 β_i z_i + ε, β_i = 1, ε ~ N(0,6)  
eps <- rnorm(n, mean = 0, sd = sqrt(6))  
# z4, z5 は Z_rest の 1 列目・2 列目  
y   <- z123[,1] + z123[,2] + z123[,3] + Z_rest[,1] + Z_rest[,2] + eps  

# —— 5. 欠損パターンの付与 ——  
# 論文の logit モデルに従い、δ1,δ2,δ3 を生成  
# δ1: logit = −1 − z4 + 2 z5 − y  
# δ2: logit = −1 − z4 + 2 z51 − y  
# δ3: logit = −1 − z50 + 2 z51 − y  
logit1 <- -1 - Z_rest[,1]    + 2*Z_rest[,2]    - y  
logit2 <- -1 - Z_rest[,1]    + 2*Z_rest[,48]   - y  
logit3 <- -1 - Z_rest[,47]   + 2*Z_rest[,48]   - y  

p1 <- plogis(logit1); p2 <- plogis(logit2); p3 <- plogis(logit3)  
δ1 <- rbinom(n, 1, p1); δ2 <- rbinom(n, 1, p2); δ3 <- rbinom(n, 1, p3)  

# δ=1 の箇所を欠損値 NA に  
z123[δ1 == 1, 1] <- NA  
z123[δ2 == 1, 2] <- NA  
z123[δ3 == 1, 3] <- NA  

# —— 6. データセット結合 ——  
sim_data <- data.frame(
  y, 
  z123,
  Z_rest
)
head(sim_data)
