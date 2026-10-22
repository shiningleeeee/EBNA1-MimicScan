# EBNA1 分子模拟挖掘：完整方法与结果报告（供外部评审）

版本 1.0.0 · 生成日期 2026-10-22 · 所有数字均可用 §14 的命令重跑复现

---

## §0 一句话结论

对 EBNA1（UniProt `P03211`，641 aa）与 20,431 条人类 reviewed 蛋白做了五条互相独立的通道（复现比对、四条零模型、表位锚定候选链、结构可及性、无偏局部相似性挖掘）。**序列层面不存在能通过组成对照的显著性信号**：所有已知模拟表位都落在低复杂度区，任何保留局部组成的零模型都能复制它们。唯一站得住的产出是一张**无偏局部相似性候选表**（363 个蛋白，其中 61 条比对长度 ≥20 aa 且一致率 ≥40%），它是一份假说排序，不是检出结果——表中最好的 E 值也只有 3.09。

---

## §1 问题定义与范围

**要回答的问题**：EBNA1 是否与某些人类蛋白存在局部序列相似性，从而可能成为抗体交叉反应（分子模拟）的靶点？哪些蛋白最值得追？

**明确不做的事**：不做抗原预测（无 HLA 分型）、不做表位预测模型、不做 docking、不声称任何候选已被验证。

**最关键的约束（决定了全部方法设计）**：EBNA1 的 N 端半段是一个极端低复杂度蛋白。任何以"相似性显著"为筛子的方法都会在这一段上产生海量假阳性，而任何以"序列不相似"为否证的论证也会被这一段毁掉。**整个工作的技术核心是如何构造一个不被低复杂度欺骗的对照**，而不是比对本身。

---

## §2 数据与溯源

| 项 | 值 |
|---|---|
| 查询序列 | EBNA1_EBVB9 `P03211`，641 aa |
| 查询 sha256 | `56d0ecacb15063e79976c62b02882ad0f49f0cb83a26a5868e51311e18ab9b0a` |
| 靶库 | Homo sapiens，reviewed (Swiss-Prot)，**20,431** 条 |
| 靶库 sha256 | `bf1bc7e188e55199a2447fc20d25834db5f7f798daa818432370deb4a6b0df5e` |
| UniProt release | **2026_03**（2026-09-02） |
| 注释快照 | `data/cache/uniprot_human_annot.tsv`，sha256 `238d9f0597427376920a082885acd82fdc8ab03b9a7a4aa234b716f48a214e0c` |
| 结构来源 | AlphaFold DB（`alphafold.ebi.ac.uk/api/prediction/{ACC}`）+ PDB 7K7R（GlialCAM 表位来源） |

`data/MANIFEST.tsv` 记录每个文件的 source_url / release / 日期 / sha256 / 序列数。流水线启动时 `s00_provenance` 校验 sha256 与 release，任一不符即抛 `ProvenanceError` 并中止——**报告里的每个数字都能追到这两个哈希**。

`data-raw/`（摘自 Biomni 交付物）另有 `benchmarks_v13.tsv`，记录五条已发表模拟主张的坐标与出处，用途见 §11。

---

## §3 计算环境

Ubuntu 22.04.5 · 14 核 · 64 GB · RTX 4090 24 GB（**GPU 全程未使用**，本工作全部是 CPU 序列/统计计算）· gcc 11.4.0 · Python 3.10.12（numpy 2.2.6 / pandas 2.3.3 / scipy 1.15.3 / biopython 1.88 / pyyaml 6.0.3 / pytest 9.1.1）· MMseqs2 静态构建 commit `fcf516aba34be4569b3eb5bda2b9db9e26a8f278`。

环境限制（影响可复现性，需知悉）：该容器 **pypi.org 不可达**（返回 000），装包走清华镜像 + uv 单二进制；系统 `apt-get install` 极慢（4 分钟未开始解包），故未安装 R/g++/cmake；venv 内无 pip。这些是整个实现选择 C 而非 R/numpy 的原因（见 §10）。

---

## §4 通道一：复现参考实现的比对扫描

### 4.1 为什么先复现

参考实现（Biomni 交付的 `ebna1-mimicscan_v1.3`）是 R + C++，其数据目录为空、无法离线复核任何数字。为了确认"我批评的到底是算法本身还是它的一次运行"，先用 C 逐位移植其比对引擎，在同一份数据上重跑。

移植遵守的语义（全部从源码逐字读出）：`kmer_selfscore = sum(diag(BLOSUM62[s,s]))`、GAr 分数 `<0.75` 的 4-mer 才作种子、UNGAPPED x-drop 6、gapped 窗口 slack 20、**严格 Smith-Waterman（v 永不为负）**、`u=H[i-1][j]-gap_open`、`l=H[i][j-1]-gap_open`（**线性间隙**）、tie 顺序 d>u>l、`extend_to_min_len()` 把扩展前最优记作 core 坐标、去重顺序 `order(-gscore)` → `unique(acc,qstart,qend)` → `unique(acc,tstart,tend)`。

### 4.2 结果

| 模式 | 本实现 | 参考实现 | 差 |
|---|---:|---:|---:|
| seeded | 59,108 | 59,092 | +16 (+0.03%) |
| exhaustive (`n_iter=4`) | 80,255 | 80,213 | +42 (+0.05%) |
| 合计 | **139,363** | **139,305** | +58 (+0.04%) |

**头部逐位一致**：LORICRIN score 542、COL1A1 475、COL3A1 473，与参考实现排第一的视图完全一致。四个已发表基准与候选均在命中集中。运行时间：seeded 1.17 s、exhaustive 25.4 s。

### 4.3 残差归因（诚实说明）

58 条差异来自去重的**同分并列**：参考实现用 `order(-gscore)` + `unique()`，同分时的保留者依赖 data.table 的内部顺序；本实现用 open-addressing 哈希，顺序确定但与 data.table 不同。**未做逐位对齐，因为源顺序不可外部推断。**

### 4.4 在此过程中发现的参考实现缺陷（与本次结论无关，但对上游有意义）

| 缺陷 | 证据 |
|---|---|
| `SEED_SCORE_MIN=14` 是 **no-op** | BLOSUM62 对角最小值为 4，任何 4-mer 自比 ≥16 > 14，实际只由 GAr 占比决定 |
| `GAP_EXT=1` 是**死参数** | DP 递推里根本没有 `gap_ext` 项，是线性间隙 |
| `MIN_NID=5` 是**死配置** | 只出现在 config、manifest 与报告，任何扫描/过滤代码都没引用 |
| `config_sha256` 不可外部复算 | 用 `digest::digest(..., algo="sha256")` 默认 `serialize=TRUE`，非 R 无法复算，而文档称其为"config.R 全文哈希" |
| `results/layers/*.tsv` 是逗号分隔 | 8 个文件后缀 `.tsv` 实为 CSV，`read.delim` 会静默读成 1 列 |
| 声称交付的三张表不存在 | `layer_counts.tsv` 声明 raw_hits / all_hits_annotated / loci_merged 三层，三表均未交付 |
| 测试数硬编码 | `18_run_manifest.R` 写 `tests_total = 115`，日志实为 128 |

---

## §5 通道二：四条互相独立的显著性检验（全部为负）

这是本工作的技术主体。"负结果"在这里是可交付的结论，因此必须证明它不是某一次实现的偶然。

### 5.1 检验 A —— 全库保组成洗牌（seeded 重扫）

把 EBNA1 整条序列洗牌（保留精确氨基酸组成），重新扫描全库。B=200：

- 每次洗牌命中数：**mean 48,096，sd 3,134，CV 0.0652**
- 真实 seeded 59,108 → **仅 1.23×**
- 真实全量 139,363 → **2.90×**

**旁证**：参考实现报告其同一零模型的 "CV = 0.0008"。但我实测 CV 是 0.0652（差 80 倍），中心值也差 18%。**它的 CV 不可能是"每次洗牌命中数"的变异系数**——这是它把零模型当作辩护性装饰的证据。

### 5.2 检验 B —— 三个并列统计量（各自配自己的零模型）

单一统计量容易被批"你选的统计量不对"，故并列三个，回归到**同一批比对**时都一样：

| 统计量 | E==0 | E<B | 中位数 | 已发表阳性落在哪 |
|---|---:|---:|---:|---|
| `E_raw`（原始 score） | 108 (0.08%) | 551 (0.40%) | 3,379,152 | — |
| `E_dens`（score/长度） | **0** | **0** | 6,666,925 | — |
| `E_band`（同长度带内） | 420 (0.30%) | 953 (0.68%) | 1,477,461 | **CRYAB 271,264 / HEPACAM 158,321 / ANO2 127,305 / MEF2B 591** |

**结论**：四个已发表的真阳性在任何一个统计量下都与噪声不可区分。**"显著性"作为筛子会把全部已知真阳性扔掉**——这是后续所有方法选择的出发点。

### 5.3 检验 C —— 真实库搜索 + 局部窗口洗牌（MMseqs2）

MMseqs2 参数：`-s 7.5 --max-seqs 20000 -e 10000 --comp-bias-corr 1 -a`。

- 真实搜索：**363 条比对，最好 E = 3.088**；E<0.001/0.01/0.1/1 全为 **0**；E<10 → 2；E<100 → 12。**四个已发表基准一条都没被报出。**
- 零模型设计要点：把真实 EBNA1 与全部洗牌**放进同一个 multi-FASTA 查询库、一次搜索**，消除批次差异（801 条查询仅 6.5 s）。B=1000：

| 零模型族 | 洗牌单位 | null best-E 中位 | P(null ≤ real) | 真实优势 |
|---|---|---:|---:|---:|
| 10-aa 窗内洗牌 | 10 aa | 12.3 | **0.176** | 4.0× |
| 20-aa 窗内洗牌 | 20 aa | 16.2 | **0.170** | 5.3× |
| 40-aa 窗内洗牌 | 40 aa | 21.4 | **0.140** | 6.9× |
| 全局洗牌 | 641 aa | 447 | 0.009 | 144.8× |

**这张表是整个结论的支点**：优势随零模型保留的局部组成的尺度单调衰减。保留 40 aa 局部组成，p=0.14，"真实 EBNA1 比随机好 145 倍"这件事完全消失。**残基顺序不携带任何可检测信号。**

### 5.4 检验 D —— 短基序穷举（6–10 mer）

不依赖任何长比对，直接穷举 EBNA1 全长的 6/7/8/9/10-mer 对全库的精确与近似匹配，每种长度配自己的窗口洗牌零模型（B=50，共 12 个可检验格子）：

| 长度 | 允许错配 | 真实命中（0/1/2 错配） | 洗牌均值 | p |
|---|---|---|---:|---:|
| 6 | 1 | 153 / 10,980 | 145.3 / 11,124.1 | 0.314 / 0.549 |
| 7 | 2 | **9** / 480 / 22,080 | 5.1 / 509.6 / 22,881.6 | **0.118** / 0.588 / 0.686 |
| 8 | 2 | 0 / 76 / 3,103 | — / 61.4 / 3,142.9 | 0.176 / 0.569 |
| 9 | 2 | 0 / 3 / 146 | — / 2.5 / 161.3 | 0.412 / 0.784 |
| 10 | 3 | **0 / 1 / 16** / 689 | — / 0.2 / 17.6 / 721.5 | 0.196 / 0.529 / 0.667 |

Bonferroni 阈值 0.05/12 = 0.0042。**最小 p = 0.118，无一显著；且真实值多数低于洗牌均值（方向也不对）。**

L=10 全库只有 17 条 ≥8/10 一致的命中，全部是组成假象（PG-rich / SR-rich / 疏水-亮氨酸）：BLTP1、CERCAM、RSPH10B/B2、TACC2、ARL13A、NFKBID、UNC80、SPG7、KMT2D、CDX2、ADGRB3×2、METTL4×2。L=7 精确 7/7 只有 9 条。

---

## §6 组成解释：为什么这四张表给出同一个答案

**EBNA1 90–330 是 241 个残基、99.5851% 为 Gly 或 Ala 的块**（G 156 / A 84 / **R 329 是唯一例外**，3 种残基；222 个 20-aa 窗口完全落在块内，窗内 GA 占比 0.95–1.00；块外最大窗口 GA 占比 0.95）。

这一段的后果是双重的：

1. **比对层**：它使参考实现的排序完全失真——其排名第一的视图 100/100 命中全部 `gar_artifact=TRUE`、全部 tier C，前 20 名是胶原（COL1A1/A2/2A1/3A1/5A1/5A2/5A3/7A1/10A1/11A1/11A2/22A1/27A1）、LORICRIN、HRNR、ELN、KRT9 等人类低复杂度蛋白。
2. **统计层**：它使任何"打乱序列再比对"的零模型失效，因为打乱后这一段仍然是最好的比对源。

**因此正确的零模型必须保留局部组成**，这就是 §5.3 那张表的全部意义。

---

## §7 通道三：表位锚定的候选链，以及它自己的零模型（自我否证）

在 §5.2 得出"显著性筛子会把真阳性全扔掉"之后，我改走另一条路：**以已发表的、被抗体/血清验证过的短表位为种子**，问"这个确切的肽段出现在哪些人类蛋白上"。

- 种子 C1 `PPRRPP`（EBNA1 394–399，d=0）、C2 `RRPFF`（402–406，d=0）、C3 `PGAIEQGPAD`（431–440，d≤3）
- 短种子**故意不放错配**：d≤1 时 C1 从 6 个蛋白暴涨到 282、C2 从 3 到 289
- 打分轴：精确匹配（2/1 分）+ CNS 表达（1）+ 膜可及（1），文献不计分

**产出 14 条**：CNTN2 5 分；MAGEE1 4；CRYAB / MEF2B / ANO2 / DOCK10 3；PRRC2A / ZNF335 / SH3PXD2B / INTS1 / MORN1 2；KCNH8 / TGOLN2 / POLE 1。

### 7.1 这个排序的零模型（这是最重要的自我否证）

把种子换成**同一残基集合的其他排列**，重跑整套打分：

| 种子 | 真实最佳 | N_A（自排列）| N_B（全库随机低复杂度肽）|
|---|---|---|---|
| C1 `PPRRPP` | 4 分 (CNTN2) | n=14，**12 个也拿到 4 分**，p=**0.8667** | p=0.1493 |
| C2 `RRPFF` | 4 分 (MAGEE1) | n=29，7 个 ≥4，p=**0.2667** | p=0.2759 |
| C3 `PGAIEQGPAD` | 3 分 (ANO2) | n=200，68 个 ≥3，p=**0.3433** | p=0.3503 |

`PPRRPP` 的 14 个排列里 12 个拿满 4 分（如 `PPPPRR`→SOCS7、`PPPRPR`→NRG2、`PPPRRP`→DDN、`PPRPPR`→ADAMTSL4）。

**结论：这套 5 分制排序没有判别力。"MAGEE1 得 4 分"描述的是人类蛋白库的成分偏倚，不是 EBNA1 的性质。** 这 14 条名单因此**不能**作为发现交付。

### 7.2 更根本的方法错误（用户指出，已更正）

这条链以 EBNA1 386–440 的三个已知表位为种子，等于**用答案锁定了搜索空间**——"照靶子画箭"。同时它与自己的论证矛盾：一边说"召回已知阳性不能评价流水线质量"，一边让已知阳性决定了名单的形状。§9 的无偏挖掘是对这一错误的直接更正。

---

## §8 通道四：结构可及性（AlphaFold 代理）

若肽段相同但埋在折叠核里，它不可能是血清学模拟。用 AlphaFold 模型算 pLDDT 与相对 SASA（Shrake-Rupley，probe 1.4 Å，92 点/原子）。

- **EBNA1 本身**：global pLDDT **41.37**，**82% 残基 pLDDT < 50**。已验证的三个表位区（394–399 / 402–406 / 431–440）全部 pLDDT 26–30 + rSASA 0.77–0.90，即**无序且高度暴露**。EBNA1 唯一 pLDDT ≥70 且长度 ≥15 的三段是 476–491（84.8）、505–539（82.9）、559–585（82.7）——DNA 结合/二聚化核心，无已知表位。
- **14 条候选**：11 条 `disordered_or_loop_exposed`，2 条 `partly_buried`（DOCK10@1550 84.1/0.24、INTS1@841 89.1/0.21），**只有 1 条 `surface_of_folded_domain`（CNTN2@911，85.1/0.26）**。
- **校准**：随机的匹配段里 `folded+exposed` 占 **19/76 = 25%**；真实候选里 1/14 = 7.1%；**p = 0.2597**。→ "CNTN2 表面可及"是常态，不是发现。

---

## §9 通道五：无偏局部相似性挖掘（当前交付物）

**设计**：query = EBNA1 全长 641 aa，不掩蔽任何一段，不以任何已知表位为种子。哪个蛋白出现，只由"EBNA1 的某一段和它像不像"决定（MMseqs2，`--comp-bias-corr 1`）。

### 9.1 结果

**363 个人类蛋白报出比对，全部 E > 1**（最好 E=3.09）。按库搜索常规标准，没有一个是显著相似的。

两套视图（答不同问题）：

- **按 E 值**（同源性搜索会报什么）：PDE11A `Q9HCR9` E=3.09 / 163 aa / 29.4% / EBNA1 464–619；MAP1S `Q66K74` E=7.08 / 46 aa / 36.9%；PLEKHH3 `Q7Z736` E=37.2 / 86 aa / 27.9%；**MAPK8IP2 `Q13387` E=64.6 / 37 aa / 54.0% / EBNA1 394–429**。
- **按局部一致率**（抗体能看到什么）：前 15 名全是 8–14 aa、70–91%，且集中在 EBNA1 541–560 的 poly-Pro/Gly 段——**组成假象，已在报告中明确标注不作发现**。

### 9.2 工作集

要求 **比对长度 ≥20 aa 且一致率 ≥40%** → **61 条**；其中落在 CNS 表达或膜可及蛋白上的 **26 条**。

| 一致率 | 长度 | E | EBNA1 | 人类蛋白 | 定位 |
|---:|---:|---:|---|---|---|
| 59.0% | 22 aa | 9.3e3 | 383–404 | COBL cordon-bleu | SPR-rich |
| 55.0% | 20 aa | 339 | 533–552 | GPC2 glypican-2（GPI 锚定） | folded core |
| 54.0% | 37 aa | 64.6 | 394–429 | MAPK8IP2 / JIP-2 | SPR-rich |
| 50.0% | 26 aa | 257 | 530–555 | FBLIM1 migfilin | folded core |
| 48.2% | 29 aa | 4.1e3 | 52–78 | RPH3A rabphilin-3A | N 端 |
| 47.6% | 21 aa | 112 | 540–560 | TSEN54 | folded core |
| 46.4% | 28 aa | 148 | 534–558 | KNDC1 | folded core |
| 45.4% | 22 aa | 777 | 535–556 | TMEM163 | folded core |
| 45.0% | 20 aa | 2.4e3 | 576–595 | PRKN parkin | folded core |
| 43.7% | 32 aa | 1.4e3 | 527–558 | SYNE1 nesprin-1 | folded core |
| 42.8% | 42 aa | 777 | 506–546 | DPH2 | folded core |

### 9.3 比对落在 EBNA1 的哪里

| 区段 | 比对数 |
|---|---:|
| folded core 476–585 | 276 |
| C-terminal | 68 |
| SPR-rich 370–445 | 12 |
| N-terminal | 7 |
| **Gly/Ala 块 90–330** | **0** |

GA 块为零不是被掩蔽，而是 MMseqs2 的组成偏差校正自己剔除了它。**这本身是一个结论**：换成无偏、组成校正的搜索，那段低复杂度块根本不构成相似性信号。

### 9.4 最值得追的一条（及其边界）

**MAPK8IP2 / JIP-2**：37 aa、54% 一致、EBNA1 394–429 ↔ JIP2 309–341，是唯一同时满足"长度够长 + 落在 SPR-rich 真实区 + 结构上折叠且表面暴露（pLDDT 85.1 / rSASA 0.26，全部候选里唯一）"的比对。JIP2 是神经元 JNK 通路支架蛋白。

**但必须同时说清**：E=64.6 在这个库规模下就是随机水平；共享肽段对交叉反应是必要不充分条件——已发表的 GlialCAM 那条是靠 7K7R 晶体结构才定下来的。

---

## §10 工程结构：可复现性

### 10.1 为什么是 C 而不是 R/numpy

三条原因：①容器里只有 gcc 11.4.0 + make，无 R、无 g++、无 cmake，apt 不可用；②需要逐位复现参考实现 C++ 引擎的 SW 语义（含 x-drop 不定长提前退出、线性间隙、tie 顺序、min-len 补长），numpy 无法向量化，纯 Python 是分钟到小时级；③只有内循环用 C，零模型/统计/绘图全 Python。

**代价（如实记录）**：700 行手写 C 与手动内存管理，调试成本高（本项目踩过 3 个 C 层面 bug：数组维度不足、静态指针未分配导致 segfault、编码序号当字符打印导致行数虚高）；去重顺序与参考实现有 58 条不可消除的残差。**若重来，会先跑 MMseqs2 拿到结论再决定是否移植。**

### 10.2 流水线

10 个 stage，单一参数源 `config.yaml`（sha256 `780bd96a27319cac26a29051c2f7c7d8af179e0e6e93854224ebe281f67ce4ed`）：

```
s00_provenance → s10_seeds → s20_motif_matches → s30_annotate → s40_priority
→ s50_alignment → s60_composition → s65_structure → s70_confirm → s80_report
```

- 任一 stage 失败即中止并写 manifest
- `reports/run_manifest.json` 记录 config 哈希 + 每个输入的 sha256 + 每个输出的 sha256 + 每阶段耗时
- s80 只从已落盘的表格插值渲染报告，不做任何计算

### 10.3 可复现性实测

| 项 | 值 |
|---|---|
| 全流程耗时 | **143.77 s**（结构模型缓存暖） |
| stage 数 | 10 |
| 输出文件 | 18 |
| 连跑两次对比 | **除 `run_manifest.json` 的 `duration_s` 外，18 个输出逐字节相同** |
| 测试 | **37 passed in 1.04 s** |
| 首次运行额外成本 | 下载约 75 MB AlphaFold 模型，约 4 min |

自洽性测试曾真实抓到一次事故：改了 `config.yaml` 却没重跑 s00，`test_config_hash_in_summary_matches_config_file` 报 summary 里是陈旧哈希——证明必须整链重跑。

### 10.4 修的三个真实 bug（记录教训）

1. `scan.py` 的 `mismatch[inside]` 索引错位 → `IndexError: boolean index did not match indexed array along axis 0; size of axis is 11438663 but size of corresponding boolean axis is 6`。正解：先 `mm = mismatch[idx]` 再按 `inside` 过滤。
2. 非确定性：`af_pdb` 无重试，一次网络抖动静默丢零模型样本 → 两次运行 p 从 0.25 变 0.26（n_null 75→76）。加三次重试，取不到的模型记入 `summary["excluded"]` 而不是当作"不可及"。
3. 取整：`summary.json` 存 0.996，测试从 FASTA 重算 0.995851，差 1.49e-4 → 全部改 6 位。

---

## §11 已发表基准的敏感性审计（不是可靠性指标）

| 蛋白 | 坐标 | 出处 | 本流水线是否"可达" |
|---|---|---|---|
| CRYAB `P02511` | 11–16 ↔ 402–407 | Thomas 2023 Sci Adv | 是（RRPFF core 5/5） |
| ANO2 `Q9NQ90` | 140–149 ↔ 431–440 | Tengvall 2019 / Lum 2026 Table 1 | 是（7/10） |
| MEF2B `Q02080` | 272–277 ↔ 394–399 | Lum 2026 PNAS（PPRRPP 小节，**非 Table 1**） | 是 |
| GlialCAM `Q14CZ8` | 370–389 ↔ 386–405 | Lanz 2022 Nature / 7K7R | **否**：只有 6/20 一致（`SSPPRA` vs `PPRRPP`），表位由晶体结构定义 |
| MBP `P02686` | 411–426 | Lum 2026 Table 1 | **否**：源文献明确是**构象型**，肽段筛法看不见 |

**召回 3/5。** 两个漏掉的都不是实现缺陷，而是方法边界。

### 11.1 关于 MEF2B 的更正（供评审注意）

参考实现把 MEF2B 登记为"held-out 阳性基准"。核查 PMID 41790925 / PMC12994170（Lum et al., PNAS 2026, 123(11):e2519445123）全文：**摘要只列三个线性模拟（GlialCAM、CRYAB、ANO2）**，Table 1 四行为 GliaCAM（线性）、CRYAB（线性）、MBP（构象）、ANO2（线性）。MEF2B 出现在二级小节 "Anti-EBNA1 antibodies also bind other proteins sharing the same linear PPRRPP motif"，是"抗 GlialCAM 克隆抗体也结合 MEF2b，因两者共享 PPRRPP"。该文同时明写 **PPRRPP 并非 EBV 独有（HHV6、部分痘病毒也有）**。

→ MEF2B 是真实的交叉反应观察，但**不是三个经典 mimic 之一**；把它升格为独立阳性基准属于过度解释（且它自述 held-out，实与其余基准同源，不独立）。

**而 PPRRPP 正好落在 EBNA1 376–404 的 SPR-rich 低复杂度区**——本文作者自己举的这个例证，同样落在组成驱动区。

---

## §12 我自己的错误清单（供评审判断可信度）

1. **在候选链上"照靶子画箭"**：以已知表位为种子锁定搜索空间。由用户指出后改为 §9 的无偏挖掘。
2. **E_band 有长度偏差**，且我把 `identity≥0.35` 任意排在显著性之后，导致"唯一幸存 CBX8"是排序选择造成的假象，不是发现。已在报告中更正为"同时通过组成过滤与校准显著性的候选，一个也没有"。
3. **曾把 EBNA1 90–330 写成"100% Gly/Ala"**，实际是 99.5851%，唯一例外是 R329。原句是 20-aa 窗口的结论被写成了逐残基陈述。
4. **组成过滤器会误杀真阳性**：HEPACAM 因 `ent_t 2.62 < 2.8` 被滤掉、MEF2B 因 `gar_q 0.63 > 0.5` 被滤掉。这是"过滤器有效"的假象来源之一。
5. **§7 的 5 分制排序**在零模型下无判别力，已自我否证。
6. 一个 12-aa 的 `PPPPRR` 之类排列就能拿满分——说明该打分轴测的是蛋白库，不是病毒。
7. 早期对参考实现 "CV=0.0008" 的批评措辞过强；正确表述是"其 CV 不可能是每次洗牌命中数的变异系数"。

---

## §13 边界与局限

1. **无一条显著**：所有通道的 p 值都在 0.12–0.87；无偏搜索最好的 E 值是 3.09。任何"候选"都是排序假说。
2. **测不到构象型表位**：MBP 型的模拟本质上需要结构/界面计算（docking、界面比对），本工作完全未做。
3. **测不到表位加工**：抗原呈递、MHC 结合、翻译后修饰（如 Lanz 2022 报道的 EBNA1 pSer376）都不在范围内。
4. **AlphaFold 模型是代理不是测量**：pLDDT 是模型置信度，rSASA 算在预测单体上，无构象系综、无翻译后修饰、无结合伴侣。
5. **库里只有 reviewed Swiss-Prot**：20,431 条。TrEMBL（约 30 倍）未查，可能漏掉未审核蛋白。
6. **"可达/不可达"的判定是基于 UniProt 文本注释的正则匹配**，会有噪音（POLE 的 121 篇 PubMed 共现就是查询污染的例子）。
7. **短基序稀有 ≠ 特异性**：`PGAIEQGPAD` 在 20,431 条人类蛋白里出现 0 次，但这不意味着它是特异性的表位。
8. **未与参考实现的全部数字逐位对齐**（58/139,363 = 0.04% 残差，原因已说明）。

---

## §14 复现命令

```bash
cd /workspace/PH-Mimic

# 完整流水线（10 stage，约 2.5 分钟，结构模型已缓存）
bash run.sh

# 单阶段 / 区间
bash run.sh --list
bash run.sh --only s65_structure

# 测试
env/venv/bin/python -m pytest tests/ -q        # 37 passed

# 无偏局部相似性挖掘（§9）
env/venv/bin/python src/mine_local_similarity.py

# 参考实现的比对复现（§4）
env/bin/mimic_scan --query data/cache/P03211.fasta \
  --db data/cache/human_swissprot.fasta --mode seeded --out reports/50_seeded.tsv
env/bin/mimic_scan --query data/cache/P03211.fasta \
  --db data/cache/human_swissprot.fasta --mode exhaustive \
  --n-iter 4 --threads 8 --out reports/50_exhaustive.tsv

# MMseqs2 真实搜索（§5.3）
env/bin/mmseqs easy-search data/cache/P03211.fasta data/cache/human_swissprot.fasta \
  reports/mmseqs_ebna1.m8 work/tmp -s 7.5 --max-seqs 20000 -e 10000 \
  --comp-bias-corr 1 -a --threads 14 --format-output \
  "target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend,qaln,taln"
```

**注意**：`--comp-bias-corr 1` 不可省（否则 GA 块会淹没结果）；`-a` 不可省（否则 convertalis 报 `Backtrace cigar is missing`）。

---

## §15 文件清单

**流水线产物（`reports/`）**

| 文件 | 内容 |
|---|---|
| `REPORT.md` | 流水线自动渲染的主报告（10 节） |
| `REPORT_local_similarity.md` | §9 的无偏局部相似性报告 |
| `run_manifest.json` | config 哈希 + 输入/输出 sha256 + 各阶段耗时 |
| `summary.json` | 全部汇总数字（自洽性测试的对象） |
| `45_local_similarity.tsv` | 363 条无偏比对原表 |
| `45_candidate_table.tsv` | 61 条 ≥20aa 且 ≥40% 的子集（含序列） |
| `50_seeded.tsv` / `50_exhaustive.tsv` | 139,363 条比对全记录（19 MB） |
| `40_null_calibration.json` / `40_seed_scores.tsv` | §7 的打分与零模型 |
| `65_structure_real.tsv` / `65_structure_null.tsv` | §8 的结构测量与校准 |
| `60_composition_profile.tsv` | 20-aa 滑动窗组成剖面 |
| `00_provenance.json` | 输入哈希与工具版本 |

**代码**：`src/phmimic/`（包：core / seeds / scan / annotate / score / structure / runner / stages/s00…s80）、`src/c`（`mimic_scan.c`、`motif_mine.c`）、`src/mine_local_similarity.py`、`src/legacy/`（step1–step5 的探索脚本）、`tests/`（37 项）。

**关键中间数据**：`reports/mmseqs_ebna1.m8`（363 条）、`reports/batch_null_per_query.csv`（B=1000 窗口洗牌）、`reports/gated_counts.tsv`（短基序 12 格零模型）、`reports/candidates_rule_R.csv`（规则 R 的无特异性验证）。

---

## §16 给评审的具体问题

1. **零模型的选择是否正确？** §5.3 显示"保留 40 aa 局部组成"后优势归零。这是否足以否证序列层面的 mimicry 信号？还是说 20–40 aa 的窗口洗牌本身就把信号洗掉了（即零模型过强），导致假阴性？如果过强，正确的替代是什么（例如保留二肽组成、保留二级结构倾向）？
2. **E<1 是否是合适的门槛？** 对表位模拟这种"短肽 + 严格"的场景，EMBOSS water 类的局部比对分数或 identity×length 是否比 BLAST 式 E 值更合适？§9 的两套视图是否应该合并成单一评分？
3. **§9 的 61 条工作集筛选（≥20 aa 且 ≥40%）是否是合理的事前规则**，还是又一次"看着答案定阈值"？
4. **结构通道在 p=0.26 的情况下是否还有信息价值**？还是说它目前只是把"候选多在无序区"这个已知事实量化了一遍？
5. **§11 把召回率定位为"敏感性审计而非可靠性指标"是否恰当**？如果评审认为应当有独立的工程可靠性判据，建议什么？
6. **是否有本文未考虑、但成本可接受的第六通道**（例如二肽/三肽组成对照、IEDB 实测表位交叉检索、TrEMBL 全库）应当先做？
7. **诚实性问题**：§12 列了我自己发现的 7 个错误。请判断其中是否有任何一条使 §0 的结论失效，或使某条通道的结果不可用。
