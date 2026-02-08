
library(data.table)
devtools::load_all()

# 1. 模拟一个包含不同品种/国家来源的系谱
# 假设:
# - A, B 是 "Breed1" (如本地品种)
# - C, D 是 "Breed2" (如引进品种)
# - E 是 A x B (纯 Breed1)
# - F 是 C x D (纯 Breed2)
# - G 是 E x F (50% Breed1, 50% Breed2)
# - H 是 G x A (75% Breed1, 25% Breed2)

ped_data <- data.table(
  Ind = c("A", "B", "C", "D", "E", "F", "G", "H"),
  Sire = c(NA, NA, NA, NA, "A", "C", "E", "G"),
  Dam = c(NA, NA, NA, NA, "B", "D", "F", "B"),
  Sex = c("male", "female", "male", "female", "male", "female", "male", "female"),
  # 只有基础个体有明确的 Breed 标签，后代的 Breed 是混合的
  Breed = c("Breed1", "Breed1", "Breed2", "Breed2", NA, NA, NA, NA)
)

tped <- tidyped(ped_data)

# 2. 计算基因流 (Ancestry Proportions)
# 这一步计算每个个体中来自 Breed1 和 Breed2 的基因比例
gf <- pedancestry(tped, labelvar = "Breed")

print("Gene Flow (Ancestry Proportions) Results:")
print(gf)

# 3. 验证结果
# G 应该是 0.5/0.5
# H 应该是 0.75/0.25 (G的0.25 + A的0.5 来自 Breed1, G的0.25 来自 Breed2)
