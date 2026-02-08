
library(visPedigree)
library(data.table)

# --- 1. 模拟一个深度系谱数据 (Simulate Pedigree) ---
# 结构：
# G0: Founders (A, B)
# G1: C (A x B), D (A x B) -> 全同胞
# G2: E (C x D) -> 全同胞交配，产生 F=0.25 的后代
# G3: F (E x A) -> 回交
# G4: G (F x E) -> 高度近交

ped_data <- data.table(
  Ind = c("A", "B", "C", "D", "E", "F", "G"),
  Sire = c(NA, NA, "A", "A", "C", "E", "E"),
  Dam = c(NA, NA, "B", "B", "D", "B", "F"),
  Sex = c("male", "female", "male", "female", "male", "female", "male"),
  Gen = c(0, 0, 1, 1, 2, 3, 4)
)

# 转换 tidyped 对象
tped <- tidyped(ped_data)

# --- 2. 基础计算 ---
# 先算一下总近交系数，作为参照
f_res <- inbreed(tped)
print("个体的总近交系数 (Total Inbreeding F):")
print(f_res[, .(Ind, f)])

# --- 3. 测试 Partial Inbreeding ---
# 我们重点关注祖先 A 和 B 对后代近交的贡献
# E (C x D): 近交来源是 A (0.125) 和 B (0.125)，总 F=0.25
# 让我们验证这一点

ancestors_to_test <- c("A", "B", "C", "D")

cat("\n--- 计算特定祖先 (A, B, C, D) 的部分近交系数 (Partial Inbreeding) ---\n")
pf_res <- pedpartial(tped, ancestors = ancestors_to_test)
print(pf_res)

# --- 4. 验证核心逻辑 ---
# 验证 E 的 pF 之和是否等于其总 F
e_pf <- pf_res[Ind == "E"]
e_total_pf <- e_pf$A + e_pf$B # C 和 D 是 Parent，但不是 Founder (Common Ancestor)
e_real_f <- f_res[Ind == "E", f]

cat(sprintf("\n个体 E 的校验:\n"))
cat(sprintf("来自 A 的 pF: %.4f\n", e_pf$A))
cat(sprintf("来自 B 的 pF: %.4f\n", e_pf$B))
cat(sprintf("pF 总和: %.4f\n", e_total_pf))
cat(sprintf("实际总 F: %.4f\n", e_real_f))

if (abs(e_total_pf - e_real_f) < 1e-6) {
  cat("✅ 校验成功：E 的部分近交之和等于总近交。\n")
} else {
  cat("❌ 校验失败\n")
}

# --- 5. 自动测试 Top 祖先功能 ---
cat("\n--- 测试自动识别 Top 祖先模式 ---\n")
# 自动识别贡献最大的前 2 个祖先
pf_auto <- pedpartial(tped, top = 2)
print(pf_auto)
