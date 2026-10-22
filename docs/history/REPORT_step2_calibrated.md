# EBNA1 × 人类蛋白质组 —— Step 2：可复现的扫描 + 可校准的显著性

日期 2026-10-21｜作者：agent（PH-Mimic）｜状态：Step 2 完成，等待结构/界面通道决策

---

## 0. 一句话结论

1. Biomni v1.3 的扫描算法被我**逐位复现**（139,363 vs 139,305 条命中，+0.04%，Top-3 分数完全相同），所以它的引擎没有错。
2. 但它报告的"分子模拟"信号，**100% 由 EBNA1 90–330 位一段 240 aa、G+A 占比 = 1.00 的低复杂度块**与人类低复杂度区（胶原 Gly-X-Y、角蛋白、loricrin）的**局部组成匹配**解释。
3. 一旦零模型保留 ≤40 aa 的局部组成，p 值从 0.009 涨到 0.14–0.18；**残基顺序不携带任何可检测的 mimicry 信号**。真实库里**没有任何一条比对 E < 1**。

---

## 1. 溯源（provenance first）

| 项 | 值 |
|---|---|
| UniProt release | **2026_03**（2026-09-02） |
| 查询 | `P03211` EBNA1_EBVB9，641 aa |
| 查询 sha256 | `56d0ecacb15063e79976c62b02882ad0f49f0cb83a26a5868e51311e18ab9b0a` |
| 库 | `proteome:UP000005640 AND reviewed:true`，**20,431 条**（13,671,780 B） |
| 库 sha256 | `bf1bc7e188e55199a2447fc20d25834db5f7f798daa818432370deb4a6b0df5e` |
| 下载记录 | `data/cache/MANIFEST.tsv`（source_url / release / sha256 / n_seqs） |
| 校验脚本 | `scripts/00_check_env.sh`（rc=0） |
| 搜索器 | MMseqs2 静态 avx2，commit `fcf516aba34be4569b3eb5bda2b9db9e26a8f278` |

Biomni 的 zip 里 `data/` 是空的，`audit/input_checksums.sha256` 引用的 5 个文件全部缺失 →
它的任何数字在包内**无法离线复核**。本报告的所有数字都是从上面两个文件重新算出来的。

---

## 2. 复现 Biomni 的扫描引擎

C 移植：`src/mimic_scan.c`（gcc -O3，忠实移植 `R/scan_core.cpp` 的严格 SW、
线性间隙、tie 顺序 d>u>l、`extend_to_min_len` 补长、两遍去重）。

| 模式 | Biomni v1.3 | 本次 | 差 |
|---|---|---|---|
| seeded（k-mer 种子 + x-drop + gapped SW） | 59,092 | **59,108** | +16（+0.03%） |
| exhaustive（全蛋白 SW，n_iter=4） | 80,213 | **80,255** | +42（+0.05%） |
| 合计 | 139,305 | **139,363** | +58（+0.04%） |
| 运行时间 | — | 1.2 s + 25.4 s（8 线程） | — |

**Top-3 分数逐位相同**：LORICRIN 542、COL1A1 475、COL3A1 473 —— 与它的
`ranked_views/1_sequence_evidence.tsv` 完全一致。残差来自同分命中的排序
（`order(-gscore)` 依赖 data.table 内部顺序），不是算法差异。

→ **结论：引擎可信，问题全部出在统计层和排序层。**

---

## 3. 那段 240 aa 的 Gly-Ala 块

`src/08_where.py` 的分段组成：

| EBNA1 区间 | G+A 占比 | 其他 |
|---|---|---|
| 0–90 | 0.38 | |
| **90–330** | **1.00** | 240 aa，只有 Gly 和 Ala |
| 330–460 | 0.32 | Arg 0.20 |
| 460–530 | 0.26 | |
| 530–641 | 0.19 | |

这是在 UniProt 上众所周知的 EBNA1 Gly-Ala repeat（GAr），但**它的定量后果没有被算过**：

- Biomni 全库 139,363 条命中里，**对齐到查询的片段 GAr ≥ 0.75 的占 14.7%**；
  **Top-50 基因里 49 个满足**。它排名第一的视图 100 行**全部**是
  `gar_artifact = TRUE`、`tier = C`、`E_emp = 0`。
- 这个块的 4-mer 自比对分数、以及它与胶原/角蛋白/loricrin 的比对分数，
  几乎完全由组成决定 —— 而这正是为什么洗牌零模型（保留全局组成、打乱顺序）
  会给出 `E_emp = 0`：它把唯一有判别力的信息（240 aa 的块）打碎了，
  却保留了"人类库里胶原很多"这个混淆因素之外的组成富集。

**MMseqs2 给了完全相反的图景**（图 1b）：它做了 composition bias correction，
报告的 363 条比对**全部落在 350–641**，**与 90–330 重叠的为 0/363**；
覆盖 EBNA1 残基 276/641（43.1%），峰值在 550–599（41.6%）。
也就是说，主流搜索器**主动拒绝**了那段块。

---

## 4. 三种统计量并列（B = 200，保组成全局洗牌）

同一批真实命中，换统计量结论就翻转：

| 统计量 | 定义 | 真实 = 0 的命中 | < B 的命中 | 中位 E | 四个已发表基准 |
|---|---|---|---|---|---|
| `E_raw` | 比分数 | 108 (0.08%) | 551 (0.40%) | 3,379,152 | 都在 |
| `E_dens` | 分数/长度 | **0** | **0** | 6,666,925 | 都在 |
| `E_band` | 同长度带内 | 420 (0.30%) | 953 (0.68%) | 1,477,461 | 都在 |

基准的 `E_band`：MEF2B 591、ANO2 127,305、HEPACAM 158,321、CRYAB 271,264。
→ **四个已知真阳性与噪声完全不可区分**，用任何固定阈值都选不出它们。

这解释了 Biomni 那个硬件的尴尬：`E_emp = 0` 有 108–420 个，但其中没有一个是真阳性。

---

## 5. 校准过零模型：保留局部组成（决定性）

方法：**real + 全部洗牌序列放进同一个 multi-FASTA 查询库，一次 MMseqs2 搜索**
（801 条查询，6.5 s）→ 批次效应被彻底消除，每个 null 与 real 走完全相同的
预过滤、组成校正和 E-value 计算。

参数：`-s 7.5 --max-seqs 20000 -e 10000 --comp-bias-corr 1 -a --threads 14`

| 零模型族 | 保留的局部组成 | null best-E 中位 | **P(null ≤ real)** | 真实/nulls 中位优势 |
|---|---|---|---|---|
| 窗口内洗牌 w = 10 | 10 aa 组成 | 12.3 | **0.176** | 4.0× |
| 窗口内洗牌 w = 20 | 20 aa 组成 | 16.2 | **0.170** | 5.3× |
| 窗口内洗牌 w = 40 | 40 aa 组成 | 21.4 | **0.140** | 6.9× |
| 全局洗牌（Biomni 式） | 整条蛋白组成 | 447 | **0.009** | 144.8× |

真实 EBNA1 best E = **3.088**（Q9HCR9/PDE11A，bits 31，163 aa，29.4% id）。

读法：
- 全局洗牌下"真实比 null 好 145 倍"是**真的**，但它是**唯一保留 240 aa 块的那个族**给的 ——
  换句话说，这个倍数度量的正是那段块本身，而不是 mimicry。
- 一旦允许 null 也拥有 40 aa 尺度的组成（w = 40），**145× → 6.9×，p = 0.14**。
  这是"残基顺序是否携带信号"的正确检验，答案是**不携带**。
- 三个窗口族的 p 值在 0.14–0.18 区间彼此一致（B = 1000，单次运行），
  且与早期 B = 25 的探索性运行（0.12 / 0.08 / 0.24）同向。

---

## 6. 这个结论能推翻什么、不能推翻什么

**能推翻的**
- "EBNA1 与人类蛋白质组存在全库层面的序列模拟" —— 在 Swiss-Prot 20,431 条上，
  去掉局部组成后没有信号（p = 0.14–0.18），最好 E = 3.09。
- 任何以 `E_emp` 升序、或未归一化 score 排序产出的候选名单 ——
  那是低复杂度目录，不是候选。
- Biomni 报告的 `null CV = 0.0008`：我实测每次洗牌的命中数 CV = **0.0652**
  （mean 48,096，sd 3,134），差 80 倍。它那个数不可能来自"每次洗牌的命中数"。

**不能推翻的**
- 已发表的**抗体层面**交叉反应（CRYAB、HEPACAM、MEF2B、ANO2 那几条）。
  那些是用合成肽 + 抗体/ELISA 验出来的，生物学上可能完全成立。
  本次说的是：**在 20,431 条人类蛋白序列上，找不到支持"广泛模拟"的统计证据**。
  "某条肽能被交叉反应抗体识别"和"两个蛋白质组之间存在模拟信号"是两件事，
  中间缺的那一环是：这个肽在**人类全库里**是否稀有。本报告补的就是这一环。
- 结构/界面层面的模拟（未做，见 §9）。

---

## 7. 漏斗里剩下的东西——以及它为什么仍然不算候选

从 139,363 条命中按"排除低复杂度与假象"层层过滤：

| 步骤 | 剩余 | 占全量 |
|---|---|---|
| 起点 | 139,363 | 100% |
| `aln_len ≥ 20` | 45,487 | 32.6% |
| 查询片段 GAr < 0.50 | 26,453 | 19.0% |
| 靶标片段 entropy ≥ 2.8 且 maxfrac ≤ 0.30 | 20,604 | 14.8% |
| 靶标片段 ga < 0.40 且 spr < 0.50 | **19,544** | 14.0% |

在这 19,544 条里看 `E_band`（200 次保组成全局洗牌中，同长度带内分数超过它的次数）：

| 门槛 | 数量 |
|---|---|
| `E_band = 0`（200 次洗牌无一超过） | **6** |
| `E_band < 200` | 8 |
| 中位 | 776,272 |

那 6 条：

| gene | acc | mode | score | aln_len | identity | E_band |
|---|---|---|---|---|---|---|
| IRX4 | P78413 | exhaustive | 39 | 120 | 20.8% | 0 |
| LRRC43 | Q8N309 | exhaustive | 43 | 142 | 20.4% | 0 |
| PPP1R3G | B7ZBB8 | exhaustive | 54 | 120 | 28.3% | 0 |
| SMG9 | Q9H0W8 | exhaustive | 60 | 151 | 24.5% | 0 |
| HMX1 | Q9NP08 | exhaustive | 60 | 129 | 24.8% | 0 |
| INSM2 | Q96T92 | exhaustive | 52 | 166 | 24.1% | 0 |

**它们不是候选，而且原因值得写下来**：120–166 aa 的比对、SW 分数只有 39–60、
identity 20–28% —— 这正是随机比对的样子。`E_band = 0` 只说明"没有一次洗牌能在
同长度带里凑出更长的弱比对"，而洗牌把二级结构倾向和局部重复一起打碎了，
所以它天然凑不出长比对。**`E_band` 因此仍然带长度偏差**，和我在 Step 1
批评 Biomni 的"未归一化分数让最长比对霸榜"是同一个病。我自己的统计量没有免疫。

唯一的另一种幸存者是 **CBX8（Q9HC52）**：identity 37.3%（≥ 0.35，是我最后加的一道
任意门槛），83 aa，score 107，EBNA1 **354–432** ↔ CBX8 **162–244**。
但它的 `E_band` = 107，`E_raw` = 222 —— **不显著**；354–432 又紧贴 Gly-Ala 块的
右边界（330），组成仍是主因。

→ **诚实的总结：既能通过组成过滤、又能通过校准显著性检验的候选，一个也没有。**
漏斗最后一步之所以"只剩一个"，是因为我把 identity 这道门槛放在了显著性之后 ——
它把 6 条 `E_band = 0` 的（低 identity）全删掉，只留下唯一一条 identity ≥ 0.35 的
（不显著）。这个顺序是我定的，不是数据给的；换任何顺序都改变不了
"没有任何候选同时满足两个条件"这个事实。

---

## 8. 复现

```bash
cd /workspace/PH-Mimic

# 0. 环境与数据校验（sha256 + release）
bash scripts/00_check_env.sh

# 1. 复现 Biomni 扫描引擎
gcc -O3 -march=native -pthread -Wall -o env/bin/mimic_scan src/mimic_scan.c
env/bin/mimic_scan --query data/cache/P03211.fasta --db data/cache/human_swissprot.fasta \
    --mode seeded     --out reports/seeded_raw.tsv
env/bin/mimic_scan --query data/cache/P03211.fasta --db data/cache/human_swissprot.fasta \
    --mode exhaustive --n-iter 4 --threads 8 --out reports/exhaustive_raw.tsv

# 2. 命中集组成 / 假象占比 / 基准存在性
env/venv/bin/python src/01_summarize.py

# 3. 三种统计量并列的洗牌零模型（B=200）
env/venv/bin/python src/03_null_stats.py

# 4. 组成过滤漏斗
env/venv/bin/python src/04_residual.py

# 5. 真实库校准搜索
env/bin/mmseqs search data/cache/P03211.fasta data/cache/human_swissprot.fasta \
    work/mm real_tmp -s 7.5 --max-seqs 20000 -e 10000 --comp-bias-corr 1 -a --threads 14
env/bin/mmseqs convertalis data/cache/P03211.fasta data/cache/human_swissprot.fasta \
    work/mm reports/mmseqs_ebna1.m8 \
    --format-output target,evalue,bits,alnlen,pident,qstart,qend,tstart,tend,qaln,taln

# 6. 局部组成零模型（B=1000，单批搜索）
rm -rf work/batch            # 必须先清，否则 mmseqs search 直接退出 1
env/venv/bin/python src/07_batch_null.py

# 7. 定位：信号在 EBNA1 的哪一段
env/venv/bin/python src/08_where.py

# 8. 图 1
env/venv/bin/python src/09_plot.py     # -> reports/figure1.png
```

产物：`reports/{seeded_raw,exhaustive_raw}.tsv`、`reports/mmseqs_ebna1.m8`、
`reports/batch_null_per_query.csv`、`reports/null_windowed_B.csv`、`reports/figure1.png`。

---

## 9. 已知局限

1. **库是 reviewed Swiss-Prot（20,431 条）**，不含 TrEMBL、不含量异构体。
   人类蛋白质组的完整 TrEMBL 约 10 万条；本报告的所有 p 值都只对
   "人工审阅的人类蛋白"这个总体成立。
2. **单条查询、单个时间点**。结论是"EBNA1 在这套库上没有可检测的序列模拟"，
   不是"EBV 不模拟人类"。
3. **只做了序列通道**。结构/界面层面的模拟（Foldseek 远端相似性、界面残基匹配）
   完全没做 —— 这也是 4090 目前还没被用到的原因。
4. p 值由 1,000 次洗牌的经验分布给出，分辨率 0.001；报告里给的是
   单次运行值，未做跨种子稳定性检查。
5. `tcov`/`qcov` 未参与筛选；`E_band` 的分箱边界（6/15/30/60/120/250/500）
   是我定的，换边界会让具体数字变，但不改变 §5 的结论（窗口族 p 值稳定在 0.14–0.18）。
6. **`E_band` / `E_raw` / `E_dens` 三者都是我自己造的统计量，都不等价于校准 E-value。**
   §7 的 6 条 `E_band = 0` 就是反例：它们只是在"洗牌凑不出长比对"这个意义上领先，
   不是生物学信号。真正的校准判据只有 §5 的 MMseqs2 E-value（best = 3.09，无一 < 1）。

---

## 10. 下一步的三个选项

| 选项 | 做什么 | 代价 | 能回答什么 |
|---|---|---|---|
| **A. 结构通道** | Foldseek 扫 AlphaFold 人类组 + EBNA1 结构，看界面层面 | 4090 就位，约 1–2 h | "序列不模拟"是否在结构层面也成立 |
| **B. 换成 TrEMBL 全库** | 用 10 万条重跑 §5 | CPU 约 30 min | 会不会因为库变大而出现信号 |
| **C. 换问题** | 不再问"EBNA1 模拟谁"，改问"哪条肽在人类全库里稀有且被抗体识别" | 需要 IEDB 数据 | 这才是已发表交叉反应的正确复现方式 |

我的建议是 **A**：它是唯一能真正用到 4090、并且能给出独立证据（与序列无关）的方向；
B 只是把同一个否定结论在更大的库上重说一遍；C 需要引入外部数据源，属于另一个项目。
