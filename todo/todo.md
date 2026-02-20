# visPedigree 项目待办事项 (TODO)


## 3. 文档与学术产出 (Documentation & Manuscript)
- [x] **完善 manuscript.md**: 
    - [x] 补充具体的算法复杂度分析（Time Complexity）。
    - [x] 补充 `deep_ped` 和 `big_family_size_ped` 的实际运行时间基准测试数据。
    - [x] 导出高清案例图示并插入文稿。

## 5. 可视化布局专项优化 (Visualization & UI)
- [ ] **边路由优化 (Edge Routing)**: 针对复杂近交系谱，探索非线性边路径，避免线条穿过无关节点。
- [ ] **垂直扩展优化**: 针对超宽系谱（> 500英寸），实现水平自动切分分页（Tiling）导出 PDF。
- [ ] **性能提升**: 将 `Forward-Backward Averaging` 等布局计算核心迁移至 `Rcpp`，并优化 `strwidth` 的调用频率。
- [ ] **交互式浏览器**: 探索集成 `htmlwidgets` (如 `D3.js`) 以支持 Web 端缩放、路径追踪 and 动态高亮。
- [ ] **时间物理轴**: 提供可选模式，将 Y 轴映射为真实出生年份，同时保持代数分层约束。

## 6. 图论和C++在包中的应用现状与性能瓶颈分析 (Graph Theory & C++ Applications and Performance Analysis)



### 6.4 潜在性能瓶颈 (Potential Performance Bottlenecks)
- [ ] **图论操作瓶颈**: 将拓扑排序和邻域搜索移至C++实现
  - 当前在R中: `topo_sort(g)`, `neighborhood()` 大图时较慢
- [ ] **内存限制**: 
  - 密集矩阵限制: >25,000个体时禁用密集矩阵
  - 稀疏矩阵构建: 大型系谱的三元组构建仍内存密集
- [ ] **Compact模式复杂度**: `compact_ped_for_matrix` 函数的O(n²)复杂度
- [ ] **R与C++交互开销**: 
  - 频繁的数据类型转换 (R vector ↔ Armadillo matrix)
  - `query_relationship` 中的映射查找可优化

### 6.5 建议的性能优化 (Suggested Performance Improvements)
#### 短期优化 (相对简单)
- [ ] **C++化图操作**: 将拓扑排序移至C++
- [ ] **预分配优化**: 更精确的向量预分配大小估算
- [ ] **缓存改进**: 扩展全同胞检测到更复杂的模式

#### 中长期优化 (需要重构)
- [ ] **流式计算**: 对超大系谱实现流式矩阵计算
- [ ] **GPU加速**: 使用CUDA/OpenCL加速大型矩阵运算
- [ ] **更智能的compact算法**: 基于图分析的更高效压缩策略

### 6.6 基于图论的深度分析扩展 (Graph-based Deep Analysis Extensions)
- [ ] **稀疏矩阵 A 矩阵计算**:
    - [ ] 利用 `Matrix` 包构建消元算子 $(I-M)$，通过前向/后向回代极速计算 $A$ 矩阵及其逆矩阵。
    - [ ] 实现针对 100 万级系谱的线性复杂度（$O(N)$）计算逻辑，替代传统的表格法。
- [ ] **点对点亲缘关系追踪 (Sargolzaei 算法)**:
    - [ ] 利用 `igraph` 的局部性原理，实时提取两个体间的最小共同祖先子图。
    - [ ] 针对大数据集实现微秒级的特定亲缘查询，无需预计算全矩阵。
- [ ] **近交系数算法优化 (Meuwissen & Luo 1992)**:
    - [ ] 结合 `data.table` 的 `by` 分组功能，实现"按家系循环"而非"按个体循环"，在 R 层级跑赢传统的 C 循环。
- [ ] **连通分支识别**: 自动识别系谱中互不相连的繁育系，帮助用户清理孤立数据或识别群体融合。
- [ ] **瓶颈节点 (Bottleneck) 检测**: 使用中介中心度 (Betweenness) 识别对群体遗传多样性至关重要的核心祖先。
- [ ] **错误自动建议**: 当检测到循环时，利用反馈弧集 (Feedback Arc Set) 算法建议最可能录错的父本/母本记录。
- [ ] **最短/全亲缘路径搜索**: 计算任意两个体间的所有遗传路径，超越简单近交系数，提供更细致的亲缘评估。

