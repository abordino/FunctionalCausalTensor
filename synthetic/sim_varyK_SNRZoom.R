setwd("~/Documents/phd/projects/causalMatrix/code")

source("mu_L4block.R")
source("oracle_mu_L4block_pool_local.R")

options(stringsAsFactors = FALSE)

library(latex2exp)

make_orthonormal = function(n, r) {
  X = matrix(rnorm(n * r), n, r)
  qr.Q(qr(X), complete = FALSE)[, 1:r, drop = FALSE]
}

make_xy = function(N2, T2, seed = 1) {
  set.seed(seed)
  x = rnorm(N2)
  y = rnorm(T2)
  x = x / sqrt(sum(x^2))
  y = y / sqrt(sum(y^2))
  list(x = x, y = y)
}

xy_to_L = function(x, y) {
  list(
    U_L = matrix(x, length(x), 1),
    Sigma_L = matrix(1, 1, 1),
    V_L = matrix(y, length(y), 1)
  )
}

true_mu = function(Mk, N1k, T1k, x, y) {
  N = nrow(Mk)
  Tt = ncol(Mk)
  Md = Mk[(N1k + 1):N, (T1k + 1):Tt, drop = FALSE]
  as.numeric(t(x) %*% Md %*% y)
}

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
  L = xy_to_L(x, y)
  
  mu_L4block_pool(
    Y = Y_list,
    k = k,
    r = r,
    U_L = L$U_L,
    Sigma_L = L$Sigma_L,
    V_L = L$V_L,
    tau = tau,
    N1 = N1_vec,
    T1 = T1_vec
  )
}

est_oracle_local = function(Y_list, k, U, V, R_list, x, y, N1_vec, T1_vec) {
  L = xy_to_L(x, y)
  
  oracle_mu_L4block_pool_local(
    Y = list(Y_list[[k]]),
    k = 1,
    U = U,
    V = V,
    R_list = list(R_list[[k]]),
    U_L = L$U_L,
    Sigma_L = L$Sigma_L,
    V_L = L$V_L,
    N1 = N1_vec[k],
    T1 = T1_vec[k]
  )
}

run_one_config = function(N, Tt, K, r, SNR_target, nrep,
                          svals = NULL, N1_range, T1_range,
                          scenario_label = "scenario",
                          base_seed = 2026, target_arm = 1,
                          tau = 2) {
  methods = c("real_pool", "oracle_local")
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

plot_mse_vs_K_pretty = function(sumdf, raw_df, file_png, main_title = NULL,
                                error_bar = c("ci95", "se", "sd", "none")) {
  error_bar = match.arg(error_bar)
  
  methods = c("oracle_local", "real_pool")
  snr_vals = sort(unique(sumdf$SNR_target))
  
  col_map = c(
    oracle_local = "#f95919",
    real_pool = "#1b9e77"
  )
  
  pch_map = c(
    oracle_local = 20,
    real_pool = 16
  )
  
  lab_map = c(
    oracle_local = TeX(r"(Oracle local $A_3$)"),
    real_pool = TeX(r"(Estimated pooled)")
  )
  
  lty_candidates = c(1, 2, 3, 4, 5, 6)
  lty_map = setNames(lty_candidates[seq_along(snr_vals)], as.character(snr_vals))
  
  Ks = sort(unique(sumdf$K))
  
  err_col = switch(
    error_bar,
    ci95 = "CI95",
    se = "SE",
    sd = "SD",
    none = NA_character_
  )
  
  if (error_bar == "none") {
    sumdf$ERR = 0
  } else {
    sumdf$ERR = sumdf[[err_col]]
    sumdf$ERR[!is.finite(sumdf$ERR)] = 0
  }
  
  pos_vals = c(
    sumdf$MSE,
    sumdf$MSE + sumdf$ERR
  )
  
  pos_vals = pos_vals[is.finite(pos_vals) & pos_vals > 0]
  eps = if (length(pos_vals) > 0) min(pos_vals) * 0.5 else 1e-16
  
  y_upper_all = log10(pmax(sumdf$MSE + sumdf$ERR, eps))
  y_lower_all = log10(pmax(sumdf$MSE - sumdf$ERR, eps))
  y_line_all = log10(pmax(sumdf$MSE, eps))
  
  y_min = min(c(y_lower_all, y_line_all), na.rm = TRUE)
  y_max = max(c(y_upper_all, y_line_all), na.rm = TRUE)
  y_rng = diff(range(y_min, y_max))
  
  if (!is.finite(y_rng) || y_rng == 0) {
    y_rng = 1
  }
  
  cu_max = max(raw_df$c_u, na.rm = TRUE)
  cell_min = min(raw_df$c_ell, na.rm = TRUE)
  
  if (is.null(main_title)) {
    main_title = TeX(r"($\log_{10}(MSE)$ vs.\ $K$)")
  } else {
    main_title = TeX(main_title)
  }
  
  error_label = switch(
    error_bar,
    ci95 = "95% Monte Carlo CI",
    se = "Monte Carlo SE",
    sd = "SD of squared errors",
    none = "No error bars"
  )
  
  png(file_png, width = 1450, height = 650, res = 140, bg = "white")
  on.exit(dev.off())
  
  par(mar = c(5, 5, 4, 15), xpd = NA, bg = "white")
  
  plot(
    Ks,
    rep(NA, length(Ks)),
    ylim = c(
      y_min - 0.08 * y_rng,
      y_max + 0.18 * y_rng
    ),
    xlab = TeX(r"($K$)"),
    ylab = TeX(r"($\log_{10}(MSE)$)"),
    main = main_title,
    type = "n"
  )
  
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "white", border = NA)
  grid(col = "white", lwd = 1.2)
  box()
  
  for (m in methods) {
    for (snr in snr_vals) {
      dm = sumdf[sumdf$method == m & sumdf$SNR_target == snr, ]
      dm = dm[order(dm$K), ]
      
      this_lty = lty_map[as.character(snr)]
      y_line = log10(pmax(dm$MSE, eps))
      
      lines(
        dm$K,
        y_line,
        type = "b",
        col = col_map[m],
        pch = pch_map[m],
        lty = this_lty,
        lwd = 2
      )
      
      if (error_bar != "none") {
        ok = is.finite(dm$ERR) & is.finite(dm$MSE) & dm$MSE > 0
        
        if (any(ok)) {
          y_low = log10(pmax(dm$MSE[ok] - dm$ERR[ok], eps))
          y_high = log10(pmax(dm$MSE[ok] + dm$ERR[ok], eps))
          
          pos = is.finite(y_low) & is.finite(y_high) & y_high > y_low
          
          if (any(pos)) {
            arrows(
              x0 = dm$K[ok][pos],
              y0 = y_low[pos],
              x1 = dm$K[ok][pos],
              y1 = y_high[pos],
              angle = 90,
              code = 3,
              length = 0.03,
              col = adjustcolor(col_map[m], alpha.f = 0.95),
              lwd = 0.8,
              lty = this_lty
            )
          }
        }
      }
    }
  }
  
  text(
    x = mean(range(Ks)),
    y = y_max + 0.13 * y_rng,
    labels = TeX(sprintf(
      "$\\max c_u = %.4f$ \\ $\\min c_{\\ell} = %.4f$ \\ Error bars: %s",
      cu_max, cell_min, error_label
    )),
    adj = c(0.5, 1),
    cex = 0.80
  )
  
  x_leg = usr[2] + 0.08 * diff(usr[1:2])
  
  legend(
    x = x_leg,
    y = usr[4],
    legend = unname(lab_map[methods]),
    col = col_map[methods],
    pch = pch_map[methods],
    lty = 1,
    lwd = 2,
    bty = "n",
    title = TeX(r"(Method)")
  )
  
  legend(
    x = x_leg,
    y = usr[4] - 0.32 * diff(usr[3:4]),
    legend = paste0("SNR = ", snr_vals),
    col = "black",
    lty = unname(lty_map[as.character(snr_vals)]),
    lwd = 2,
    bty = "n",
    title = "SNR"
  )
}

N = 100
Tt = 80
r = 6
Ks = c(50, 150, 300, 500, 1000)

SNR_targets = c(1, 100, 10000)
svals = seq(2, 0.6, length.out = r)

nrep = 500
base_seed = 231198
tau = 0.01

N1_window = c(70, 70)
T1_window = c(60, 60)
scenario_label = "window"

raw_all = vector("list", length(SNR_targets))
idx_all = 1

for (snr in SNR_targets) {
  cat(sprintf("\nRunning SNR_target = %s\n", snr))
  
  raw_list = vector("list", length(Ks))
  
  for (i in seq_along(Ks)) {
    K = Ks[i]
    print(K)
    
    raw_list[[i]] = run_one_config(
      N = N,
      Tt = Tt,
      K = K,
      r = r,
      SNR_target = snr,
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
  
  raw_all[[idx_all]] = do.call(rbind, raw_list)
  idx_all = idx_all + 1
}

raw_df = do.call(rbind, raw_all)
sumdf = summarize_mse(raw_df)

cat(sprintf("\nOverall max c_u   = %.6f\n", max(raw_df$c_u, na.rm = TRUE)))
cat(sprintf("Overall min c_ell = %.6f\n", min(raw_df$c_ell, na.rm = TRUE)))

for (snr in SNR_targets) {
  tmp = raw_df[raw_df$SNR_target == snr, ]
  
  cat(sprintf(
    "SNR_target=%s: max c_u = %.6f, min c_ell = %.6f\n",
    snr,
    max(tmp$c_u, na.rm = TRUE),
    min(tmp$c_ell, na.rm = TRUE)
  ))
}

print(sumdf)

dir.create("synthetic/result", recursive = TRUE, showWarnings = FALSE)
dir.create("synthetic/figure", recursive = TRUE, showWarnings = FALSE)

snr_tag = paste(SNR_targets, collapse = "_")

tag = paste0(
  "multiSNR_", snr_tag,
  "_N1win", N1_window[1], "-", N1_window[2],
  "_T1win", T1_window[1], "-", T1_window[2]
)

raw_file = file.path(
  "synthetic/result",
  paste0("Zoomraw_compare_pooling_oracleLocal_varyK_randomXY_", tag, ".csv")
)

sum_file = file.path(
  "synthetic/result",
  paste0("Zoomsummary_compare_pooling_oracleLocal_varyK_randomXY_", tag, ".csv")
)

fig_file = file.path(
  "synthetic/figure",
  paste0("Zoomcompare_logMSE_vs_K_pooling_oracleLocal_randomXY_", tag, ".png")
)

write.csv(raw_df, raw_file, row.names = FALSE)
write.csv(sumdf, sum_file, row.names = FALSE)

plot_mse_vs_K_pretty(
  sumdf = sumdf,
  raw_df = raw_df,
  file_png = fig_file,
  main_title = sprintf(
    "$\\log_{10}(MSE)$ vs. $K$ (SNR targets = %s, $N_1=%d$, $T_1=%d$)",
    paste(SNR_targets, collapse = ", "),
    N1_window[1],
    T1_window[1]
  ),
  error_bar = "ci95"
)

## Reload saved results and regenerate plot without rerunning simulations
raw_file_2 = "synthetic/4Block/result/Zoomraw_compare_pooling_oracleLocal_varyK_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.csv"
sum_file_2 = "synthetic/4Block/result/Zoomsummary_compare_pooling_oracleLocal_varyK_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.csv"

raw_df_loaded_2 = read.csv(raw_file_2)
sumdf_loaded_2 = read.csv(sum_file_2)

print(sumdf_loaded_2)

fig_file_loaded_2 = file.path(
  "synthetic/4Block/figure",
  "Zoomcompare_logMSE_vs_K_pooling_oracleLocal_randomXY_multiSNR_1_100_10000_N1win70-70_T1win60-60.png"
)

plot_mse_vs_K_pretty(
  sumdf = sumdf_loaded_2,
  raw_df = raw_df_loaded_2,
  file_png = fig_file_loaded_2,
  main_title = sprintf(
    "$\\log_{10}(MSE)$ vs. $K$ (SNR targets = %s, $N_1=%d$, $T_1=%d$)",
    paste(SNR_targets, collapse = ", "),
    N1_window[1],
    T1_window[1]
  ),
  error_bar = "ci95"
)