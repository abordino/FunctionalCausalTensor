setwd("~/Desktop/code")

source("bilinearTensorAllFunction.R")

options(stringsAsFactors = FALSE)

make_orthonormal = function(n, r) {
  X = matrix(rnorm(n * r), n, r)
  qr.Q(qr(X), complete = FALSE)[, seq_len(r), drop = FALSE]
}

make_xy = function(N2, T2, seed) {
  set.seed(seed)
  
  x = rnorm(N2)
  y = rnorm(T2)
  
  list(
    x = x / sqrt(sum(x^2)),
    y = y / sqrt(sum(y^2))
  )
}

true_mu = function(Mk, N1k, T1k, x, y) {
  Md = Mk[
    seq.int(N1k + 1L, nrow(Mk)),
    seq.int(T1k + 1L, ncol(Mk)),
    drop = FALSE
  ]
  
  drop(crossprod(x, Md %*% y))
}

# The oracle estimator might not be well defined when the relevant Gram
# matrices are not invertible. This function return NA when this is the 
# case; later on, we will remove NAs when aggregating the results.
safe_mu = function(expr) {
  z = tryCatch(expr, error = function(e) NA_real_)
  z = as.numeric(z)
  if (!is.finite(z)) NA_real_ else z
}

compute_c_from_subspaces = function(U, V, N1_vec, T1_vec) {
  K = length(N1_vec)
  
  cU_lo = numeric(K)
  cU_hi = numeric(K)
  cV_lo = numeric(K)
  cV_hi = numeric(K)
  
  N = nrow(U)
  Tt = nrow(V)
  
  for (k in 1:K) {
    U1k = U[1:N1_vec[k], , drop = FALSE]
    V1k = V[1:T1_vec[k], , drop = FALSE]
    
    eu = eigen(crossprod(U1k), symmetric = TRUE, only.values = TRUE)$values
    ev = eigen(crossprod(V1k), symmetric = TRUE, only.values = TRUE)$values
    
    cU_lo[k] = min(eu) / (N1_vec[k] / N)
    cU_hi[k] = max(eu) / (N1_vec[k] / N)
    cV_lo[k] = min(ev) / (T1_vec[k] / Tt)
    cV_hi[k] = max(ev) / (T1_vec[k] / Tt)
  }
  
  list(
    c_ell = min(c(cU_lo, cV_lo)),
    c_u = max(c(cU_hi, cV_hi))
  )
}

gen_panel = function(N, Tt, K, r, N1_vec, T1_vec, SNR_target,
                     seed = 1, svals = NULL) {
  set.seed(seed)
  
  U = make_orthonormal(N, r)
  V = make_orthonormal(Tt, r)
  
  if (is.null(svals)) {
    svals = exp(seq(log(20), log(0.6), length.out = r))
  }
  
  svals = as.numeric(svals)
  min_smin = min(abs(svals))
  
  M_list = vector("list", K)
  R_list = vector("list", K)
  
  for (k in 1:K) {
    A = make_orthonormal(r, r)
    B = make_orthonormal(r, r)
    Rk = A %*% diag(svals, r, r) %*% t(B)
    R_list[[k]] = Rk
    M_list[[k]] = U %*% Rk %*% t(V)
  }
  
  sigma2 = as.numeric(SNR_target) * min_smin^2 / N
  sigma = sqrt(sigma2)
  
  Y_list = vector("list", K)
  
  for (k in 1:K) {
    E = matrix(rnorm(N * Tt, 0, sigma), N, Tt)
    Y = M_list[[k]] + E
    Y[(N1_vec[k] + 1):N, (T1_vec[k] + 1):Tt] = NA
    Y_list[[k]] = Y
  }
  
  list(
    Y_list = Y_list,
    M_list = M_list,
    U = U,
    V = V,
    R_list = R_list,
    sigma2 = sigma2
  )
}

make_N1T1_vec = function(N, Tt, K, N1_range, T1_range, seed = 1) {
  clip_int = function(x, lo, hi) max(lo, min(hi, as.integer(x)))
  
  set.seed(seed)
  
  N1_lo = clip_int(N1_range[1], 1, N - 1)
  N1_hi = clip_int(N1_range[2], 1, N - 1)
  T1_lo = clip_int(T1_range[1], 1, Tt - 1)
  T1_hi = clip_int(T1_range[2], 1, Tt - 1)
  
  draw = function(lo, hi, size) {
    vals = seq.int(lo, hi)
    vals[sample.int(length(vals), size = size, replace = TRUE)]
  }
  
  list(
    N1_vec = draw(N1_lo, N1_hi, K),
    T1_vec = draw(T1_lo, T1_hi, K)
  )
}

est_real_pool = function(Y_list, k, r, x, y, tau, N1_vec, T1_vec) {
  mu_xy4block_pool(
    Y = Y_list,
    k = k,
    r = r,
    x = x,
    y = y,
    tau = tau,
    N1 = N1_vec,
    T1 = T1_vec
  )
}

est_real_nopool = function(Y_list, k, r, x, y, tau, N1_vec, T1_vec) {
  mu_xy4block_pool(
    Y = list(Y_list[[k]]),
    k = 1L,
    r = r,
    x = x,
    y = y,
    tau = tau,
    N1 = N1_vec[k],
    T1 = T1_vec[k]
  )
}

est_oracle_pool = function(Y_list, k, U, V, R_list, x, y,
                           N1_vec, T1_vec) {
  oracle_mu_xy4block_pool(
    Y = Y_list,
    k = k,
    U = U,
    V = V,
    R_list = R_list,
    x = x,
    y = y,
    N1 = N1_vec,
    T1 = T1_vec
  )
}

est_oracle_nopool = function(Y_list, k, U, V, R_list, x, y,
                             N1_vec, T1_vec) {
  oracle_mu_xy4block_pool(
    Y = list(Y_list[[k]]),
    k = 1L,
    U = U,
    V = V,
    R_list = list(R_list[[k]]),
    x = x,
    y = y,
    N1 = N1_vec[k],
    T1 = T1_vec[k]
  )
}

est_oracle_local = function(Y_list, k, U, V, R_list, x, y,
                            N1_vec, T1_vec) {
  oracle_mu_xy4block_pool_local(
    Y = list(Y_list[[k]]),
    k = 1L,
    U = U,
    V = V,
    R_list = list(R_list[[k]]),
    x = x,
    y = y,
    N1 = N1_vec[k],
    T1 = T1_vec[k]
  )
}

run_one_config = function(N, Tt, K, r, SNR_target, nrep,
                          svals = NULL, N1_range, T1_range,
                          scenario_label = "scenario",
                          base_seed = 2026, target_arm = 1,
                          tau = 2) {
  methods = c(
    "real_pool",
    "real_nopool",
    "oracle_pool",
    "oracle_nopool",
    "oracle_local"
  )
  
  rows = vector("list", length(methods) * nrep)
  idx = 1
  
  for (rep_id in 1:nrep) {
    seed_rep = base_seed + 100 * rep_id + 100 * K + as.integer(round(SNR_target))
    
    NT = make_N1T1_vec(
      N = N,
      Tt = Tt,
      K = K,
      N1_range = N1_range,
      T1_range = T1_range,
      seed = seed_rep + 7777
    )
    
    N1_vec = NT$N1_vec
    T1_vec = NT$T1_vec
    
    N1k = N1_vec[target_arm]
    T1k = T1_vec[target_arm]
    N2k = N - N1k
    T2k = Tt - T1k
    
    dat = gen_panel(
      N = N,
      Tt = Tt,
      K = K,
      r = r,
      N1_vec = N1_vec,
      T1_vec = T1_vec,
      SNR_target = SNR_target,
      seed = seed_rep,
      svals = svals
    )
    
    cc = compute_c_from_subspaces(dat$U, dat$V, N1_vec, T1_vec)
    
    xy = make_xy(N2k, T2k, seed = seed_rep + 999)
    x = xy$x
    y = xy$y
    
    mu_t = true_mu(dat$M_list[[target_arm]], N1k, T1k, x, y)
    
    mu_rp = safe_mu(
      est_real_pool(
        Y_list = dat$Y_list,
        k = target_arm,
        r = r,
        x = x,
        y = y,
        tau = tau,
        N1_vec = N1_vec,
        T1_vec = T1_vec
      )
    )
    
    mu_r1 = safe_mu(
      est_real_nopool(
        Y_list = dat$Y_list,
        k = target_arm,
        r = r,
        x = x,
        y = y,
        tau = tau,
        N1_vec = N1_vec,
        T1_vec = T1_vec
      )
    )
    
    mu_op = safe_mu(
      est_oracle_pool(
        Y_list = dat$Y_list,
        k = target_arm,
        U = dat$U,
        V = dat$V,
        R_list = dat$R_list,
        x = x,
        y = y,
        N1_vec = N1_vec,
        T1_vec = T1_vec
      )
    )
    
    mu_o1 = safe_mu(
      est_oracle_nopool(
        Y_list = dat$Y_list,
        k = target_arm,
        U = dat$U,
        V = dat$V,
        R_list = dat$R_list,
        x = x,
        y = y,
        N1_vec = N1_vec,
        T1_vec = T1_vec
      )
    )
    
    mu_oL = safe_mu(
      est_oracle_local(
        Y_list = dat$Y_list,
        k = target_arm,
        U = dat$U,
        V = dat$V,
        R_list = dat$R_list,
        x = x,
        y = y,
        N1_vec = N1_vec,
        T1_vec = T1_vec
      )
    )
    
    mu_map = list(
      real_pool = mu_rp,
      real_nopool = mu_r1,
      oracle_pool = mu_op,
      oracle_nopool = mu_o1,
      oracle_local = mu_oL
    )
    
    for (m in methods) {
      mu_hat = mu_map[[m]]
      
      rows[[idx]] = data.frame(
        scenario = scenario_label,
        N = N,
        T = Tt,
        r = r,
        K = K,
        SNR_target = SNR_target,
        sigma2 = dat$sigma2,
        N1k = N1k,
        T1k = T1k,
        c_ell = cc$c_ell,
        c_u = cc$c_u,
        rep = rep_id,
        method = m,
        mu_true = mu_t,
        mu_hat = mu_hat,
        sq_err = (mu_hat - mu_t)^2
      )
      
      idx = idx + 1
    }
  }
  
  do.call(rbind, rows)
}

summarize_mse = function(df) {
  agg = aggregate(
    sq_err ~ scenario + N + T + r + K + SNR_target + method,
    data = df,
    FUN = function(z) {
      nn = sum(!is.na(z))
      
      if (nn == 0) {
        return(c(
          mean = NA_real_,
          sd = NA_real_,
          se = NA_real_,
          ci95 = NA_real_,
          n = 0
        ))
      }
      
      mse = mean(z, na.rm = TRUE)
      sdz = sd(z, na.rm = TRUE)
      sez = sdz / sqrt(nn)
      ci95z = 1.96 * sez
      
      c(
        mean = mse,
        sd = sdz,
        se = sez,
        ci95 = ci95z,
        n = nn
      )
    }
  )
  
  data.frame(
    scenario = agg$scenario,
    N = agg$N,
    T = agg$T,
    r = agg$r,
    K = agg$K,
    SNR_target = agg$SNR_target,
    method = agg$method,
    MSE = agg$sq_err[, "mean"],
    SD = agg$sq_err[, "sd"],
    SE = agg$sq_err[, "se"],
    CI95 = agg$sq_err[, "ci95"],
    n = agg$sq_err[, "n"]
  )
}

#----------------------------
# Run the simulation
#----------------------------

N = 100
Tt = 80
r = 6
Ks = c(1, 2, 5, 10, 20, 50, 200)

SNR_target = 1
svals = seq(2, 0.6, length.out = r)

nrep = 500
base_seed = 231198
tau = 0.01

N1_window = c(70, 70)
T1_window = c(60, 60)
scenario_label = "window"

raw_list = vector("list", length(Ks))

for (i in seq_along(Ks)) {
  K = Ks[i]
  print(K)
  
  raw_list[[i]] = run_one_config(
    N = N,
    Tt = Tt,
    K = K,
    r = r,
    SNR_target = SNR_target,
    nrep = nrep,
    svals = svals,
    N1_range = N1_window,
    T1_range = T1_window,
    scenario_label = scenario_label,
    base_seed = base_seed,
    target_arm = 1,
    tau = tau
  )
}

raw_df = do.call(rbind, raw_list)
summary_df = summarize_mse(raw_df)

results_dir = "Results"
dir.create(results_dir, showWarnings = FALSE)

tag = paste0(
  "varyK_randomXY",
  "_N", N,
  "_T", Tt,
  "_r", r,
  "_SNR", SNR_target,
  "_N1", N1_window[1], "-", N1_window[2],
  "_T1", T1_window[1], "-", T1_window[2],
  "_nrep", nrep,
  "_tau", tau
)

raw_file = file.path(results_dir, paste0("raw_", tag, ".csv"))
summary_file = file.path(results_dir, paste0("summary_", tag, ".csv"))

write.csv(raw_df, raw_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)

cat("Saved:\n", raw_file, "\n", summary_file, "\n")