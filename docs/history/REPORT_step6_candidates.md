# EBNA1 分子模拟候选 — 第 6 轮：撤掉一条假说，装上一道新过滤器

**日期**：2026-10-22
**上游**：`REPORT_candidates_v3_honest.md`、`PROJECT_DOSSIER.md`、GPT 第二轮复审任务单
**本轮新增**：成熟度过滤器（`src/51_processing_state.py`）、结构层（`src/52_biology_first.py`）、分层表（`src/53_layered_table.py`）

---

## §0 一句话结论

上一轮排名第一的 **MMP19 被撤回**——它的匹配片段落在会被切掉的前肽里。同时发现，我上一轮用来把 6 条压成 1 条的"结构可及性"轴**本身是无效的**：它会把四个已发表的 EBNA1 表位全部排除。本轮交付 22 条候选，按序列排序轴（`margin_z`）排列，各生物学轴只作标注、不再相加。

---

## §1 MMP19 的撤回

GPT 的 P0-2 指出 MMP19 是需加工的前体，我据此查了 UniProt Q99542（entry version 228，2026-09-02）：

| 特征 | 位置 | 说明 |
|---|---|---|
| SIGNAL | 1..18 | 信号肽 |
| **PROPEP** | **19..97** | **前肽——43–66 完整落在这里** |
| CHAIN | 98..508 | 成熟酶 |
| MOTIF | 83..90 | 半胱氨酸开关 |
| PTM | — | Activated by autolytic cleavage after Lys-97 |

匹配片段 `EGSNNFKPEDITEALRAFQEASEL`（前体 43–66）在激活时随前肽一起被切掉。**成熟 MMP19 上不存在这 24 个残基。**

GPT 要求区分的三个命题，答案分别是：

| 命题 | 结论 |
|---|---|
| 前体模型 rSASA 高、折叠暴露 | **成立**（AlphaFold 建模的是含信号肽+前肽的全长前体，pLDDT 73.3） |
| 切割后该片段仍存在 | **不成立**（成熟酶上不存在） |
| 天然蛋白中抗体可及 | **未知**，需要证明 proMMP19（而非成熟酶）在可及表面以足够丰度存在 |

补充一条我不知道的张力：UniProt 的 TISSUE SPECIFICITY 写 "Not detected in brain and peripheral blood leukocytes"，与我在上一版引用的"激活 PBMC 表面表达"（PMID 9562866，1998）和"活动性 MS 病灶高表达"（PMID 17083473，2006）需要分辨——正常组织 vs 激活/病灶状态可以并存，但我没有读这两篇原文，不应把它们说成"三条独立支持"。

**处理**：MMP19 从候选名单移除，作为"曾提出的假说 + 它为什么不成立"记录在案。

---

## §2 新增：成熟度过滤器

思路：序列相似只是候选的必要条件之一；匹配片段必须在抗体能遇到的蛋白形式里存在。UniProt 标注的信号肽、前肽、转运肽都在成熟过程中被切掉。

数据：`data/cache/uniprot_processing.tsv`（20,431 条 human reviewed，字段含 Signal peptide / Propeptide / Chain / Transit peptide / Initiator methionine / Peptide / Subcellular location / Tissue specificity），一次 stream 下载 13 秒、12 MB。

363 条候选的分类结果：

| 分类 | 条数 |
|---|---|
| retained_in_chain | 350 |
| retained_unannotated | 1 |
| **spans_removed_boundary** | **4** |
| **fully_in_propep** | **3**（含 MMP19） |
| **fully_in_signal** | **3** |
| **fully_in_transit** | **2** |

**12/363 条被这道过滤拦下。**

**这道过滤的校准**：三个已验证阳性在人类蛋白侧的片段——CRYAB 11–16、ANO2 140–149、GlialCAM 370–389——**全部 retained**（0/3 假阴性）。这是本轮唯一一道既非平凡、又不在已知阳性上出错的过滤器。

---

## §3 被证伪的轴：EBNA1 侧的"折叠且暴露"

上一轮我写道"结构轴把 6 条压成 1 条，只有 MMP19 两侧都折叠且暴露"。这个推理有两个问题。

**第一个问题**（GPT 已提）：预测结构可及 ≠ 抗体可识别。

**第二个问题**（本轮发现，更致命）：**用"EBNA1 片段必须折叠"来筛选，会把四个已发表的表位全部排除。**

| 已发表 EBNA1 表位 | pLDDT | rSASA | `classify` 结果 |
|---|---|---|---|
| 386–405（GlialCAM） | 28.2 | 0.83 | disordered_or_loop_exposed |
| 394–399（MEF2B/PPRRPP） | 30.4 | 0.90 | disordered_or_loop_exposed |
| 402–406（CRYAB/RRPFF） | 27.0 | 0.77 | disordered_or_loop_exposed |
| 431–440（ANO2） | 26.0 | 0.81 | disordered_or_loop_exposed |

四个全在无序区（pLDDT 26–30），全部 < 70 的折叠阈值。

**所以 "EBNA1 侧折叠且暴露" 这道轴的召回是 0/4。** 它不能当过滤器，也不能当排序轴。上一轮它之所以"把 6 条压成 1 条"，正是因为 EBNA1 的折叠区（476–585，即 DNA 结合/二聚化域）与已验证表位所在的 N 端无序区完全不重叠——这条轴把候选**推向了已知表位不会出现的地方**。

这一条推翻的是我上一轮的核心论证，不只是 MMP19 这一条。

---

## §4 分区可达性：只能标注，不能过滤

另一个自然的轴：靶蛋白必须在抗体能到达的区室（分泌/质膜胞外）。55 条里有 12 条满足。

**但它会排除 CRYAB——CRYAB 是胞质内小热激蛋白，却是三个已验证 mimic 之一。** 所以这道轴有已知假阴性，只能作为标注。

顺带说明：这条也解释了为什么不能把"胞内 ⇒ 不可能"写进结论。MS 病灶里有组织破坏与蛋白释放，胞内蛋白成为抗体靶点是有先例的（CRYAB 本身就是）。

---

## §5 结果：22 条候选（`margin_z ≥ 6`）

过滤链：mmseqs2 无偏全长搜索 363 条 → 长度 ≥20 aa 且 identity ≥40% → 成熟度 retained → 逐对比对保组成零模型 `margin_z ≥ 6`。**四轴分别列出，不相加。**

| # | 基因 | 靶标 | z | 长度 | id% | CNS | SURF | 可及区室 | EBNA1 侧 | 靶标侧 | EBNA1 区间 | 靶标区间 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **PRKAR1B** | P31321 | **10.56** | 24 | 41.6 | – | Y | Y | partly_buried | surface_folded | 564–587 | 20–43 |
| 2 | ZNF347 | Q96SE7 | 9.14 | 21 | 42.8 | – | – | – | partly_buried | disordered | 566–586 | 230–250 |
| 3 | RNF228 | A0A7I2V3R4 | 9.10 | 21 | 52.3 | – | – | – | disordered | disordered | 63–83 | 207–227 |
| 4 | CARD17P | Q5XLA6 | 8.95 | 25 | 44.0 | – | – | – | surface_folded | partly_buried | 578–601 | 56–80 |
| 5 | R3HDM4 | Q96D70 | 8.58 | 26 | 46.1 | – | – | – | surface_folded | surface_folded | 478–503 | 169–193 |
| 6 | RFLNB | Q8N5W9 | 8.51 | 27 | 40.7 | – | – | – | disordered | disordered | 588–614 | 78–104 |
| 7 | ACACA | Q13085 | 8.35 | 24 | 41.6 | Y | – | – | surface_folded | surface_folded | 488–511 | 1938–1958 |
| 8 | DALRD3 | Q5D0E6 | 8.34 | 21 | 42.8 | – | – | – | disordered | surface_folded | 401–421 | 424–444 |
| 9 | BLOC1S1 | P78537 | 8.12 | 24 | 45.8 | – | – | – | partly_buried | disordered | 516–539 | 10–33 |
| 10 | SYNE1 | Q8NF91 | 7.72 | 32 | 43.7 | Y | – | – | disordered | not_modelled | 527–558 | 5267–5297 |
| 11 | HECTD3 | Q5T447 | 7.69 | 23 | 43.4 | – | – | – | disordered | partly_buried | 584–606 | 342–363 |
| 12 | TSEN54 | Q7Z6J9 | 7.28 | 21 | 47.6 | Y | – | – | disordered | disordered | 540–560 | 257–277 |
| 13 | RHBDF2 | Q6PJF5 | 7.25 | 22 | 45.4 | – | Y | Y | disordered | disordered | 541–562 | 339–360 |
| 14 | ACACB | O00763 | 7.16 | 27 | 40.7 | – | – | – | surface_folded | surface_folded | 488–514 | 2049–2072 |
| 15 | **MAPK8IP2** | Q13387 | 6.90 | **37** | **54.0** | Y | Y | Y | disordered | disordered | 394–429 | 309–341 |
| 16 | IRF3 | Q14653 | 6.47 | 35 | 40.0 | – | – | – | disordered | partly_buried | 580–606 | 360–393 |
| 17 | PCDH8 | O95206 | 6.39 | 24 | 41.6 | Y | Y | Y | disordered | disordered | 530–553 | 868–891 |
| 18 | CCHCR1 | Q8TD31 | 6.33 | 22 | **54.5** | – | – | – | surface_folded | surface_folded | 482–502 | 143–164 |
| 19 | PRKN | O60260 | 6.31 | 20 | 45.0 | Y | – | – | surface_folded | surface_folded | 576–595 | 408–425 |
| 20 | FBLIM1 | Q8WUP2 | 6.30 | 26 | 50.0 | Y | – | – | disordered | disordered | 530–555 | 135–160 |
| 21 | WDFY3 | Q8IZQ1 | 6.24 | 27 | 40.7 | Y | – | – | partly_buried | not_modelled | 557–583 | 311–336 |
| 22 | TMEM179 | Q6ZVK1 | 7.53 | 21 | 42.8 | – | – | – | surface_folded | surface_folded | 477–497 | 72–92 |

（`margin_z = (obs_score - null_mean) / null_sd`，零模型 = 两侧片段各自保组成洗牌，B=10⁴，与 `mimic_scan.c` 同一套 SW。）

---

## §6 三条值得细看的

### 6.1 PRKAR1B — 序列轴第一，且有一个具体的结构故事

```
EBNA1  564-587  VFLQTHIFAEVLKDAIKDLVMTKP
PRKAR1B 20-43   LYVQLHGIQQVLKDCIVHLCISKP
                ..V.L.H..Q.V.L.KD.I..H.L..I.KP     10/24 identical
```

两侧都是**二聚化螺旋**：EBNA1 564–587 在其 C 端二聚化/DNA 结合域；PRKAR1B 20–43 在其 N 端 dimerization/docking 域。相同的不是散在的残基，而是 **V-L-K-D-I-L 的疏水七肽节律**。

支持：`margin_z` 10.56 是 363 条里最高；UniProt 标注 PRKAR1B 定位 `Cell membrane`（PubMed 23115245、27028580）；UniProt 的 `FT REGION 2..136 /note="Dimerization and phosphorylation"` 确认 20–43 落在二聚化区（该区还含 `REGION 67..98 Disordered`，说明二聚化区本身部分是柔性的）；PRKAR1B 在脑内表达，其突变导致 PSP 样神经退行（这是疾病相关性，不是自身免疫证据）；自身免疫语境文献计数 1。

必须写明的反面：我的 `ax_cns` 没有点亮它（UniProt 的组织特异性文本里没有脑相关关键词），所以"脑表达"这一条来自我的领域知识而非本表字段；PRKAR1B 的经典定位是胞质，`Cell membrane` 那两条文献我没有读；两侧都是螺旋，而螺旋的疏水七肽节律本身就会产生中等程度的成对一致率，所以"结构故事"不能替代显著性。

### 6.2 MAPK8IP2 — 唯一落在已验证表位区间上的

```
EBNA1   394-429  PPRRPPPGRRPFFHPVGEADYFEYHQEGGPDGEPDV
MAPK8IP2 309-341 PPREPPRRPAFLPVGPDDTNSEYESGSESEPDL
```

- 它覆盖了 EBNA1 394–429，**同时包含 PPRRPP（394–399）与 RRPFF（402–406）两个已验证表位**；
- 37 aa、54.0% identity，是 22 条里最长且 identity 最高的第二条；
- 靶标 JIP2 是脑/胰岛表达的 JNK 支架蛋白，UniProt 注明 "accumulates in cell surface projections"。

反面：MMseqs2 给它的 E=64.6，在这个库规模就是随机水平；两侧都是无序区；`PPREPPRRPA` 里的 PPR 重复很可能是脯氨酸富集的组成巧合。**它是"结构与位置都对得上但绝对相似度不显著"的典型。**

### 6.3 PCDH8 — 两侧都无序但属于黏附分子家族

原钙粘蛋白 8，单次跨膜 I 型膜蛋白，树突/突触前。EBNA1 530–553 vs 868–891，24 aa / 41.6%。胞外域是大段无序，所以"两侧无序"在这里不奇怪，也不等于不可及。

---

## §7 各轴的能力边界（写给自己和评审）

| 轴 | 能做什么 | 已知假阴性 | 已知假阳性 |
|---|---|---|---|
| `margin_z` 保组成零模型 | 排除"组成能解释"的匹配 | — | 不归一化长度，跨长度不可严格比较 |
| 成熟度 retained | **唯一 3/3 召回的过滤器** | 未知 | 无注释的蛋白默认通过（1 条 retained_unannotated） |
| 分区可及 | 标注优先级 | **CRYAB（胞质）** | 注释关键词粗粒度 |
| EBNA1 侧折叠+暴露 | **不能用**——4/4 假阴性 | 全部已知表位 | — |
| 靶标侧折叠+暴露 | 弱标注 | CRYAB 11–16 在柔性 N 端 | 预测结构 ≠ 溶液构象 |

---

## §8 相对上一版的变化

| 项 | 上一版 | 本版 |
|---|---|---|
| 头部候选 | MMP19（8 分） | PRKAR1B（z=10.56） |
| MMP19 | 第一 | **撤回**（前肽） |
| 结构轴 | 用来把 6 条压成 1 条 | **判定无效**（4/4 假阴性） |
| 过滤链 | 8 轴加权求和 | 四轴分列，不相加 |
| 成熟度 | 无 | **新增，12/363 被拦，校准 3/3** |
| 名单规模 | 6 | 22（z≥6），完整 55 |

---

## §9 复现

```bash
# 一次下载带成熟度特征的 UniProt 注释（13 s / 12 MB）
curl -s -G https://rest.uniprot.org/uniprotkb/stream \
  --data-urlencode "query=organism_id:9606 AND reviewed:true" \
  --data-urlencode "format=tsv" \
  --data-urlencode "fields=accession,gene_names,length,ft_signal,ft_propep,ft_chain,ft_transit,ft_init_met,ft_peptide,cc_subcellular_location,cc_tissue_specificity,keyword" \
  -o data/cache/uniprot_processing.tsv

env/venv/bin/python src/51_processing_state.py    # -> reports/51_processing_state.tsv
env/venv/bin/python src/52_biology_first.py       # -> reports/52_biology_first.tsv（需 AlphaFold 模型，冷缓存约 6 min）
env/venv/bin/python src/53_layered_table.py       # -> reports/53_layered_table.tsv
```

版本锚点：UniProt 2026_03、human reviewed 20,431 条、AlphaFold DB（2026-10-22 缓存）、MMseqs2 `fcf516ab…`。

---

## §10 我认为下一步该做什么

1. **PRKAR1B 值得单独查**：它是本轮唯一一个"序列信号最强 + 有具体结构解释（螺旋七肽节律）+ 有疾病相关性"的候选，而且它从未在之前的 8 轴排名里出现过——因为旧打分被 CNS 关键词和组成轴埋掉了。应该核实 PRKAR1B 的膜定位文献（23115245、27028580），并查它是否在 MS 病灶或脑脊液中被检测到。
2. **MAPK8IP2 需要一次能证伪的实验设计**：它覆盖已验证表位，若成立是最直接的"扩展"；若不成立，说明 PPRRPP 区段的所有匹配都是组成巧合。
3. **不要再加轴了。** 现在的轴已经互相矛盾（CRYAB 会同时被两道轴排除）。真正能推进的是外部证据：IEDB 上有没有针对这 22 个蛋白的实测抗体记录。
