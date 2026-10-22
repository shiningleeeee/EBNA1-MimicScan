# Step 4 — 候选名单（特异性优先筛选）

生成时间：2026-10-22 · 远端 `/workspace/PH-Mimic`
数据：UniProt release 2026_03（2026-09-02），EBNA1 `P03211`（641 aa），Homo sapiens reviewed 20,431 条

## 0. 一句话

**之前所有的筛法都从"全库相似度"出发，被组成假象淹没。这一步反过来：从抗体实验真正验证过的 EBNA1 表位出发，问"哪些人类蛋白带着这个精确基序"，再用生物学（中枢神经表达 / 膜表面可及 / 已知自身抗原）排序——这一维任何组成统计量都给不出来。结果 14 个蛋白，其中 4 个是已发表的，10 个是新的。**

## 1. 筛法

**种子（EBNA1 坐标为 1-based 闭区间，全部落在 386–440 这段被抗体验证过的交叉反应区内）：**

| 种子 | 序列 | EBNA1 坐标 | 来源 | 错配预算 |
|---|---|---|---|---|
| C1 | `PPRRPP` | 394–399 | Lum 2026 PNAS（GlialCAM 单抗交叉结合 MEF2B） | **0** |
| C2 | `RRPFF` | 402–406 | Thomas 2023 Sci Adv（CRYAB，5/5 核心） | **0** |
| C3 | `PGAIEQGPAD` | 431–440 | Tengvall 2019 PNAS（ANO2，7/10） | 3 |

**为什么短种子不开放错配**：把 C1 放宽到 1 个错配，人类组里命中数从 6 个蛋白暴涨到 282 个（`PPPRPP`、`PPRLPP`、`PTRRPP` 这类脯氨酸串）；C2 从 3 个涨到 289 个。那就是 Step 2/Step 3 已经证明会产出假象的同一个区间。所以错配预算只花在够长的种子上（C3 是 10-mer，且已发表的 ANO2 本身就是 7/10，不给 3 个错配会把它自己漏掉）。

**排序轴**（按序）：
1. 精确匹配数（`n_ident == core_len`）
2. 生物学支持分（中枢神经表达 / 膜或分泌可及 / 有疾病关联，来自 UniProt `SUBCELLULAR LOCATION`、`TISSUE SPECIFICITY`、`DISEASE`、关键词）
3. 一致率

## 2. 名单

### A 档 — 精确基序 + 有独立生物学支持

| # | 蛋白 | 基因 | 基序证据 | UniProt 定位 / 表达 | 独立证据 |
|---|---|---|---|---|---|
| 1 | Contactin-2 / TAG-1 `Q02246` | CNTN2 | `PPRPP`… 精确 `PPRRPP` @911 (6/6) | Cell membrane；中枢神经轴突 | **已知 MS 自身抗原**（Derfuss 2009 PNAS, PMID 19416878：contactin-2 定向自身免疫在 MS 患者中检出，介导灰质病理） |
| 2 | Alpha-crystallin B `P02511` | CRYAB | 精确 `RRPFF` @11 (5/5) | 胞质/核；晶状体及其他组织 | **已发表 mimic**（Thomas 2023 Sci Adv, OR 2.0 / 合并 9.0） |
| 3 | MEF2B `Q02080` | MEF2B | 精确 `PPRRPP` @272 (6/6) | 核；骨骼肌/心肌 | **已发表**（Lum 2026 PNAS：抗 GlialCAM 单抗交叉结合 MEF2b） |
| 4 | Anoctamin-2 `Q9NQ90` | ANO2 | `PGDIELGPLD` @140 (7/10) | Cell membrane；视网膜/神经 | **已发表 mimic**（Tengvall 2019 PNAS, OR 1.6） |
| 5 | MAGE-E1 `Q9HCI5` | MAGEE1 | 精确 `RRPFF` @710 (5/5) | Cell membrane / 核周 | 癌-睾丸抗原；同基序还出现在 264、348 位 |

### B 档 — 精确基序，生物学支持尚待验证

| # | 蛋白 | 基因 | 基序证据 | 备注 |
|---|---|---|---|---|
| 6 | PRRC2A / BAT2 `P48634` | PRRC2A | 精确 `PPRRPP` @1706 (6/6) | **位于 MHC class III 区**（MHC 与 MS 强关联） |
| 7 | Zinc finger protein 335 `Q9H4Z2` | ZNF335 | 精确 `RRPFF` @720 (5/5) | 核内转录辅因子；神经发育相关（小头畸形基因） |
| 8 | SH3PXD2B `A1X283` | SH3PXD2B | 精确 `PPRRPP` @758 (6/6) | 足体/细胞突起；成纤维细胞 |
| 9 | Integrator complex subunit 1 `Q8N201` | INTS1 | 精确 `PPRRPP` @841 (6/6) | 核；snRNA 加工 |
| 10 | MORN1 `Q5T089` | MORN1 | 精确 `PPRRPP` @18 (6/6) | N 端，功能注释少 |

### C 档 — 仅 7/10，基序信息量低（`Gx-EQGP` 骨架）

| # | 蛋白 | 基因 | 基序 | 共同骨架 | 备注 |
|---|---|---|---|---|---|
| 11 | DOCK10 `Q96BY6` | DOCK10 | `PSAFFQGPAD` @1550 | G..EQGP.. | 树突棘；脑/肺低表达 |
| 12 | KCNH8 `Q96L42` | KCNH8 | `TGAHEQNPAD` @973 | G..EQGP.. | 电压门控钾通道；主要在神经系统 |
| 13 | TGOLN2 / TGN46 `O43493` | TGOLN2 | `SGAEEQGPID` @224 | G..EQGP.. | 反式高尔基网络膜蛋白 |
| 14 | POLE `Q07864` | POLE | `PGAIRDGPAT` @1304 | G..EQGP.. | DNA 聚合酶 ε |

C 档的共同骨架是 `Gx-xEQGP`，这 4 条更像"ANO2 那个 7/10 的可容忍替代"，不像独立信号。**建议只把 A/B 档（10 个）当作候选，C 档留作敏感性记录。**

## 3. 这个筛法为什么不是又一次组成假象

- 它**不**先在全库排序再取 top（Step 2 的错误）；它先固定一个实验验证过的短表位，再问"谁带着它"。
- 精确 6-mer/5-mer 在 20,431 条蛋白里只命中 6 个和 3 个——集合小到可以逐个读，而不是 7,072 个。
- 排序的第二维（CNS 表达、膜可及、已知自身抗原）**与氨基酸组成完全独立**。这是它能突破组成天花板的原因。

## 4. 诚实交代：这个筛法的召回与漏检

- **召回 3/4 已发表线性 mimic**：CRYAB、MEF2B（经 C1）、ANO2（经 C3）。
- **漏掉 GlialCAM `Q14CZ8`**：它和 EBNA1 386–405 的 20-mer 只有 **6/20（30%）** 一致，没有 `PPRRPP`（它有 `SSPPRA`）。也就是说 Lanz 2022 那个用晶体结构（7K7R）定义的 mimic **不是线性可捕获的**——这正是本文自己把它标成"composition-level"的原因。
- 因此：**这个筛法测不到构象型/组成型 mimic，只测得到线性精确基序。** 拿它当"完整性判据"是错的，拿它当"线性候选发生器"是可以的。
- 6-mer 本身的绝对稀有度不等于特异性：随机挑一个 PR-rich 的 6-mer 也能得到几个蛋白。所以 A 档的排序**主要靠生物学那一维**，不是靠基序分数。

## 5. 文件

| 文件 | 内容 |
|---|---|
| `src/12_candidate_cores.py` | 种子扫描（含 d≤1 的对照，证明放宽即崩） |
| `src/13_shortlist.py` | 精确种子 + UniProt 注释 + 分档 |
| `reports/candidates_cores.tsv` | 567 个蛋白（d≤1 版本，作为"放宽后崩坏"的记录） |
| `reports/candidates_shortlist.tsv` | 14 个候选 + 注释列 |
| `reports/candidates_cores_detail.json` | 前 30 条逐基序明细 |

## 6. 建议的下一步（按性价比）

1. **对 A 档 5 条做逐条文献核查**：每条查"该蛋白 + autoantibody / multiple sclerosis / EBV"是否已有报道。CNTN2 已知；CRYAB/ANO2/MEF2B 已知；MAGEE1 待查。
2. **对 B 档 5 条查表达与可及性**：是否在 CNS 表达（HPA/GTEx）、是否膜外可及。若都不满足，降级。
3. **不要**再回去调全库相似度的阈值——三条独立零模型已经一致否证那条路。
