# tests/manual/test_timevar_formats.R
# 验证 pedgenint/pedstats 的不同 timevar 列格式与 cycle_length 参数
# 运行方式：在项目根目录下执行 devtools::load_all(".") 后 source 此文件

library(data.table)

# ── 构造基础小系谱 ───────────────────────────────────────────────
# G0 (4 founders) -> G1 (4 individuals) -> G2 (3 individuals)
base_ped <- data.frame(
  Ind  = c("G0_S1","G0_S2","G0_D1","G0_D2",
           "G1_A","G1_B","G1_C","G1_D",
           "G2_A","G2_B","G2_C"),
  Sire = c(NA, NA, NA, NA,
           "G0_S1","G0_S1","G0_S2","G0_S2",
           "G1_A","G1_B","G1_B"),
  Dam  = c(NA, NA, NA, NA,
           "G0_D1","G0_D2","G0_D1","G0_D2",
           "G1_C","G1_C","G1_D"),
  Sex  = c("male","male","female","female",
           "male","male","female","female",
           "male","female","male"),
  stringsAsFactors = FALSE
)

# ── Case 1: 纯整数年份 integer year (最常见) ─────────────────────
# G0=2018, G1=2020, G2=2022 => 世代间隔均为 2 年
cat("=== Case 1: Integer Year, unit='year' ===\n")
ped1 <- base_ped
ped1$BirthYear <- c(2018,2018,2018,2018,
                    2020,2020,2020,2020,
                    2022,2022,2022)
tped1 <- suppressMessages(tidyped(ped1))
gi1 <- suppressMessages(pedgenint(tped1, timevar = "BirthYear", unit = "year"))
print(gi1)
stopifnot(
  all(c("Pathway","N","Mean","SD") %in% names(gi1)),
  abs(gi1[Pathway == "Average", Mean] - 2.0) < 1e-9
)
cat("PASS: Average Mean == 2.0 years\n\n")

# ── Case 2: ISO 日期字符串 "YYYY-MM-DD"，以天为单位 ──────────────
# G0=2020-01-01, G1=2020-06-15 (165d later), G2=2020-12-01 (169d after G1)
cat("=== Case 2: ISO date string 'YYYY-MM-DD', unit='day' ===\n")
ped2 <- base_ped
ped2$HatchDate <- c("2020-01-01","2020-01-01","2020-01-01","2020-01-01",
                    "2020-06-15","2020-06-15","2020-06-15","2020-06-15",
                    "2020-12-01","2020-12-01","2020-12-01")
tped2 <- suppressMessages(tidyped(ped2))
gi2 <- suppressMessages(pedgenint(tped2, timevar = "HatchDate", unit = "day"))
print(gi2)
stopifnot(
  attr(gi2, "unit") == "day",
  !"GenEquiv" %in% names(gi2)  # 未指定 cycle_length 时不应有此列
)
cat(sprintf("Average Mean = %.1f days\n", gi2[Pathway == "Average", Mean]))
cat("PASS: unit='day', no GenEquiv column\n\n")

# ── Case 3: 非标准日期格式 "DD/MM/YYYY"，须指定 format= ─────────
# 与 Case 2 同一批数据，只是格式不同
cat("=== Case 3: Custom date format 'DD/MM/YYYY', unit='day', format='%d/%m/%Y' ===\n")
ped3 <- base_ped
ped3$HatchDate <- c("01/01/2020","01/01/2020","01/01/2020","01/01/2020",
                    "15/06/2020","15/06/2020","15/06/2020","15/06/2020",
                    "01/12/2020","01/12/2020","01/12/2020")
tped3 <- suppressMessages(tidyped(ped3))
gi3 <- suppressMessages(pedgenint(tped3, timevar = "HatchDate", unit = "day",
                                   format = "%d/%m/%Y"))
print(gi3)
stopifnot(
  isTRUE(all.equal(gi2$Mean, gi3$Mean, check.attributes = FALSE))
)
cat("PASS: Case 2 (ISO) and Case 3 (custom format) produce identical Mean values\n\n")

# ── Case 4: cycle_length 参数 —— 世代等效数 GenEquiv ─────────────
# 对虾等水产物种：目标世代周期为 165 天
# GenEquiv = Mean / cycle_length
# G0->G1 约 165 天 => GenEquiv ≈ 1.0
cat("=== Case 4: cycle_length=165, expect GenEquiv = Mean/165 ===\n")
gi4 <- suppressMessages(pedgenint(tped2, timevar = "HatchDate", unit = "day",
                                   cycle_length = 165))
print(gi4)
stopifnot(
  "GenEquiv" %in% names(gi4),
  isTRUE(all.equal(gi4$GenEquiv, gi4$Mean / 165, check.attributes = FALSE))
)
cat(sprintf("Average GenEquiv = %.3f (Mean=%.1fd / cycle=165d)\n",
            gi4[Pathway == "Average", GenEquiv],
            gi4[Pathway == "Average", Mean]))
cat("PASS: GenEquiv == Mean / cycle_length\n\n")

# ── Case 5: 以月为单位，float 年份列 ────────────────────────────
# 若系谱原始年份是 numeric（如 2020.5 = 2020年6月），unit 仍为 year 时保持原值
cat("=== Case 5: Numeric float year, unit='year' ===\n")
ped5 <- base_ped
ped5$Year <- c(2020.0, 2020.0, 2020.0, 2020.0,
               2020.5, 2020.5, 2020.5, 2020.5,
               2021.0, 2021.0, 2021.0)
tped5 <- suppressMessages(tidyped(ped5))
gi5 <- suppressMessages(pedgenint(tped5, timevar = "Year", unit = "year"))
print(gi5)
stopifnot(
  abs(gi5[Pathway == "Average", Mean] - 0.5) < 1e-9
)
cat("PASS: Average Mean == 0.5 years for 6-month intervals\n\n")

# ── Case 6: pedstats 整合测试 ────────────────────────────────────
cat("=== Case 6: pedstats with timevar='HatchDate', unit='day', cycle_length=165 ===\n")
st <- suppressMessages(pedstats(tped2, timevar = "HatchDate", unit = "day",
                                 cycle_length = 165))
stopifnot(
  !is.null(st$gen_intervals),
  "GenEquiv" %in% names(st$gen_intervals),
  attr(st$gen_intervals, "unit") == "day"
)
cat("gen_intervals:\n")
print(st$gen_intervals)
cat("PASS: pedstats correctly passes timevar/unit/cycle_length to pedgenint\n\n")

cat("========================================\n")
cat("All 6 cases PASSED.\n")
cat("========================================\n")
