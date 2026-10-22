# Step 5 — 候选名单（加入正交生物学轴后的重排）

生成时间：2026-10-22 · 远端 `/workspace/PH-Mimic`
数据：UniProt release 2026_03 · Human Protein Atlas（`proteinatlas.org/{ENSG}.json`）· NCBI E-utilities

## 0. 一句话

Step 4 只按"精确基序"排，会漏掉生物学上更值得追的东西。加上**中枢神经表达**和**膜外可及**这两个与氨基酸组成完全无关的轴之后，名单重排：**DOCK10 从 C 档（序列只有 7/10）升到与已发表 mimic 同分**，而几个"精确 6/6 但脑内不表达、胞内定位"的条目降到 B 档。

## 1. 排序规则（事先声明，5 分制）

| 维度 | 分值 | 判据 |
|---|---|---|
| 基序精确度 | 2 / 1 | 精确匹配 6/6 或 5/5 → 2；仅 7/10 → 1 |
| CNS 表达 | 1 | HPA `RNA tissue specific nTPM['brain'] ≥ 20` **且** 脑表达 cluster 属神经元/少突胶质/星形胶质/小胶质 |
| 膜外可及 | 1 | HPA Secretome = Secreted/Cell membrane，或 UniProt 定位为 Cell membrane / Secreted / Extracellular |
| 已有文献 | 1 | PubMed 有该基因 ×（autoantibody / autoantigen / multiple sclerosis / Epstein-Barr）共现 |

## 2. 最终名单（按总分降序）

| 分 | 蛋白 | 基因 | 基序证据 | CNS（HPA brain nTPM / 脑 cluster） | 可及性 | 文献 | 状态 |
|---|---|---|---|---|---|---|---|
| **5** | Contactin-2/TAG-1 `Q02246` | CNTN2 | `PPRRPP` @911 **6/6** | **157.3** / Oligodendrocytes–Myelination | Secreted in brain | 4 篇；**已知 MS 自身抗原**（Derfuss 2009 PNAS, PMID 19416878） | 已知 |
| **4** | MAGE-E1 `Q9HCI5` | MAGEE1 | `RRPFF` @710 **5/5**（另有 264、348） | **25.0** / Neurons | Cell membrane | **0 篇** | **新增** |
| 3 | αB-crystallin `P02511` | CRYAB | `RRPFF` @11 **5/5** | 低 / Schwann–Melanocyte | Secreted & membrane | 27 篇 | 已发表 |
| 3 | MEF2B `Q02080` | MEF2B | `PPRRPP` @272 **6/6** | 低 / Plasma cells | 胞内 | 5 篇 | 已发表 |
| 3 | Anoctamin-2 `Q9NQ90` | ANO2 | `PGDIELGPLD` @140 **7/10** | 低 / Endothelial | Cell membrane | 12 篇 | 已发表 |
| 3 | DOCK10 `Q96BY6` | DOCK10 | `PSAFFQGPAD` @1550 **7/10** | **34.2** / Oligodendrocytes–Myelination | 胞内 | DOCK-D 家族见于小鼠神经炎症模型（PMID 32241915） | **新增，被生物学抬升** |
| 2 | PRRC2A/BAT2 `P48634` | PRRC2A | `PPRRPP` @1706 **6/6** | 低 / Non-specific | 胞内 | 2 篇；**MHC class III 区** | **新增** |
| 2 | ZNF335 `Q9H4Z2` | ZNF335 | `RRPFF` @720 **5/5** | 低（神经发育基因，但成人脑内 mRNA 低） | 核 | 0 篇 | **新增** |
| 2 | SH3PXD2B `A1X283` | SH3PXD2B | `PPRRPP` @758 **6/6** | 低 / Non-specific | 足体 | 0 篇 | **新增** |
| 2 | INTS1 `Q8N201` | INTS1 | `PPRRPP` @841 **6/6** | 低 / Non-specific | 核 | 1 篇 | **新增** |
| 2 | MORN1 `Q5T089` | MORN1 | `PPRRPP` @18 **6/6** | 24.4 但 cluster 为纤毛细胞（非神经） | 胞内 | 0 篇 | **新增** |
| 1 | KCNH8 `Q96L42` | KCNH8 | `TGAHEQNPAD` @973 7/10 | 10.2 / 垂体 | 膜 | 0 篇 | 新增 |
| 1 | TGOLN2/TGN46 `O43493` | TGOLN2 | `SGAEEQGPID` @224 7/10 | 低 / Macrophage | 高尔基膜 | 0 篇 | 新增 |
| 1 | POLE `Q07864` | POLE | `PGAIRDGPAT` @1304 7/10 | 低 / Neutrophil | 核 | 121 篇（多为 MS 遗传学研究，非自身抗原证据） | 新增 |

**真正的新增候选（去掉 4 个已发表 + 3 个纯序列偶然）= 7 条**：
MAGEE1（最高分新增）、DOCK10、PRRC2A、ZNF335、SH3PXD2B、INTS1、MORN1

## 3. 这一步改了什么

- **DOCK10 上来了**：序列只有 7/10（`PSAFFQGPAD` vs `PGAIEQGPAD`），但它是少突胶质/髓鞘表达 + 胞内 GEF + DOCK-D 家族已在神经炎症模型中被研究。纯序列筛法会把它埋在几百条里，生物学轴把它挑出来——这正是序列统计量给不出的信息。
- **MORN1 降了**：brain nTPM 24.4 看起来达标，但 HPA 的 cluster 是"纤毛细胞"（室管膜纤毛），不是神经细胞。所以 CNS 判据要两条都满足，不能只看数字。
- **ZNF335 降了**：教科书式的神经发育基因，但成人脑内 mRNA 很低（HPA 未达阈值），核内定位也不便抗体接触。

## 4. 边界（不变）

- 这个筛法测不到构象型/组成型 mimic（GlialCAM 只有 6/20 一致，会被漏掉）。
- 6-mer 的稀有度本身不等于特异性；名单的排序**主要靠 CNS/可及性/文献**这三条与组成无关的轴。
- POLE 的 121 篇是查询污染（MS 遗传学 GWAS 里常有 POLE），**不能**当作自身抗原证据。

## 5. 文件

| 文件 | 内容 |
|---|---|
| `src/13_shortlist.py` | 精确种子筛选 + UniProt 注释 |
| `src/14_annotate_hpa.py` | HPA Ensembl 解析（注意：UniProt 给的是带版本的 `ENSGxxxx.x`，HPA 只认不带版本的形式）+ 脑表达 + PubMed 共现计数 |
| `reports/candidates_shortlist.tsv` | 14 条基础表 |
| `reports/candidates_annotated.tsv` | 加上 HPA/PubMed 注释后的表 |

## 6. 下一步

1. **MAGEE1** 值得优先验证：精确 `RRPFF`、神经元表达、膜蛋白、无任何自身免疫文献——是名单里"证据缺口最大但可得性最好"的一条。
2. **CNTN2** 已有 MS 自身抗原文献，价值在于提供一个**正对照**：如果抗体交叉反应假设成立，抗 EBNA1 386–405 抗体应当也结合 contactin-2，这个可以用现有试剂直接做。
3. 若要继续计算而非实验：对 A 档 4 条做 EBNA1 386–440 与候选片段的**结构可及性**检查（是否位于柔性连接区 vs 折叠域内部），这决定线性表位是否真能被抗体看到。
