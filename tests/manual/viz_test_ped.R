
library(visPedigree)
library(data.table)

# 复现测试中的系谱数据
ped_data <- data.table(
  Ind = c("A", "B", "C", "D", "E", "F", "G"),
  Sire = c(NA, NA, "A", "A", "C", "E", "E"),
  Dam = c(NA, NA, "B", "B", "D", "B", "F"),
  Sex = c("male", "female", "male", "female", "male", "female", "male"),
  Gen = c(0, 0, 1, 1, 2, 3, 4)
)

tped <- tidyped(ped_data)

# 绘制系谱图
# compact = TRUE 可以让图更紧凑，但在这么小的系谱里也没关系
visped(tped, compact = TRUE, show_graph = TRUE)
