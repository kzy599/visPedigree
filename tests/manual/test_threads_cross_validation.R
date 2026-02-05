
devtools::load_all()
library(data.table)

# --- 1. 数据准备 ---
# 我们需要一个适中大小的系谱 (N ~ 5000-10000)
# 太小看不出并行优势，太大 (N > 25000) 会触发 D 矩阵的内存保护限制

cat(">>> 准备测试数据 (目标 N=6000)...\n")
# 从 small_ped 开始扩充，而不是用巨大的 big_family_size_ped
base_ped <- tidyped(small_ped) # small_ped 大概几十个个体
target_n <- 6000

# 扩充系谱
ped_list <- list()
n_copies <- ceiling(target_n / nrow(base_ped))

for (i in 1:n_copies) {
  p <- copy(base_ped)
  p[, Ind := paste0(Ind, "_", i)]
  # 处理非空父母
  p[Sire != "", Sire := paste0(Sire, "_", i)]
  p[Dam != "", Dam := paste0(Dam, "_", i)]
  ped_list[[i]] <- p
}
tped <- tidyped(rbindlist(ped_list))

cat(sprintf("测试用系谱大小: %d 个体\n", nrow(tped)))
cat("注意: N < 25000 以避免触发生 D 矩阵的内存保护限制。\n")
cat("----------------------------------------------------------------\n")

# --- 2. 定义测试函数 ---

run_comparison <- function(label, method, invert_col = "auto", n_rep = 3) {
  cat(sprintf("\n=== 测试场景: %s (method='%s', invert='%s') ===\n", label, method, invert_col))
  
  run_timing <- function(thr) {
    times <- numeric(n_rep)
    cat(sprintf("  Running threads=%d: ", thr))
    for (i in 1:n_rep) {
      pt <- proc.time()
      # invisible() 防止打印矩阵
      invisible(pedmat(tped, method = method, invert_method = invert_col, threads = thr, compact = FALSE, sparse = TRUE))
      times[i] <- (proc.time() - pt)["elapsed"]
      cat(".")
    }
    mean_time <- mean(times)
    cat(sprintf(" Avg: %.4fs\n", mean_time))
    return(mean_time)
  }
  
  t1 <- run_timing(1)
  t4 <- run_timing(4)
  
  speedup <- t1 / t4
  is_affected <- if (speedup > 1.2) "YES (显著加速)" else if (speedup < 0.9) "NO (变慢)" else "NO (无变化)"
  
  cat(sprintf("  结果: Speedup = %.2fx -> 受 threads 参数影响? %s\n", speedup, is_affected))
  return(data.frame(
    Scenario = label,
    Method = method,
    Invert = invert_col,
    Time_1T = t1,
    Time_4T = t4,
    Speedup = speedup,
    Affected = is_affected
  ))
}

# --- 3. 执行交叉测试 ---

results <- list()

# CASE 1: D Matrix (预期：受控，加速)
results[[1]] <- run_comparison("D Matrix (Dominance)", "D")

# CASE 2: Ainv Matrix (预期：受控，加速 - 前提是 N > 5000)
results[[2]] <- run_comparison("Ainv Matrix (Inverse A)", "Ainv")

# CASE 3: A Matrix (预期：不受控)
results[[3]] <- run_comparison("A Matrix (Additive)", "A")

# CASE 4: f (Inbreeding) (预期：不受控)
results[[4]] <- run_comparison("Inbreeding (f)", "f")

# CASE 5: AA Matrix (预期：不受控，由BLAS决定)
results[[5]] <- run_comparison("AA Matrix (Epistatic)", "AA")

# CASE 6: Dinv (Auto/General) (预期：不受控)
# 注意：Dinv 计算非常慢，如果 N 很大可能跑很久，这里如果 N>2000 可能要小心
if (nrow(tped) > 3000) {
    cat("\n[Skip Dinv/AAinv tests because N is large and inversion O(n^3) takes too long]\n")
} else {
    results[[6]] <- run_comparison("Dinv (Invert D)", "Dinv", "auto")
}


# --- 4. 汇总报告 ---
final_df <- do.call(rbind, results)
cat("\n================================================================\n")
cat("                       测试结果汇总报告                        \n")
cat("================================================================\n")
print(final_df, row.names = FALSE)
cat("================================================================\n")
cat("结论验证:\n")
cat("1. 'D' 和 'Ainv' 应显示加速 (Speedup > 1.0)。\n")
cat("2. 'A', 'f', 'AA', 'Dinv' 应无显著变化 (Speedup ~ 1.0)。\n")
