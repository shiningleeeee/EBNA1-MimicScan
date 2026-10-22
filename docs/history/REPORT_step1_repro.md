# PH-Mimic · Step 1 — 复现 Biomni 的 EBNA1 mimic 扫描，并补上它缺的零模型

日期：2026-10-21 · 主机：106.75.244.36 `/workspace/PH-Mimic` · 数据：UniProt release 2026_03（2026-09-02）
查询：EBNA1_EBVB9 `P03211`，641 aa，sha256 `56d0ecac…9b0a`
库：Homo sapiens reviewed（Swiss-Prot）20,431 条，sha256 `bf1bc7e1…df5e`

---

## 1. 复现结果：算法层逐位一致

`src/mimic_scan.c` 是从 `R/seeds.R` + `R/scan.R` + `R/scan_core.cpp` 逐行移植的 C 实现
（远端无 R、无 g++，故用 gcc）。参数与 Biomni v1.3 `config/config.R` 完全相同：
`KMER=4, SEED_SCORE_MIN=14, SEED_GAR_MAX=0.75, XDROP=6, SLACK=20, GAP_OPEN=11, GAP_EXT=1(死参数), MIN_LEN=6, MIN_SCORE=20, EXHAUSTIVE_ITER=4`。

| 层 | Biomni v1.3 | 本复现 | 差 |
|---|---|---|---|
| seeded 命中 | 59,092 | **59,108** | +16 (0.03%) |
| exhaustive 命中 | 80,213 | **80,255** | +42 (0.05%) |
| 合计原始命中 | 139,305 | **139,363** | +58 (0.04%) |

耗时：seeded 1.2 s / exhaustive 25 s（8 线程）。

残差来源：`scan_proteome()` 的 `order(-gscore)` + 两次 `unique()` 在**同分时的先后顺序**
依赖 data.table 的内部顺序，我的 open-addressing 去重给出确定性但非逐位的等价结果。
除此之外，头部命中的分数**逐位相同**：

| 基因 | Biomni (ranked view 1) | 本复现 |
|---|---|---|
| LORICRIN | score 542 | **542** |
| COL1A1 | 475 | **475** |
| COL3A1 | 473 | **473** |

三个已发表基准也都在：CRYAB(P02511) 6 条 / HEPACAM(Q14CZ8) 9 条 / MEF2B(Q02080) 12 条 / ANO2(Q9NQM4) 12 条。

**结论：Biomni 的比对层是真实可复现的，不是报告层面的修饰。**

---

## 2. 它没交付的东西，现在有了：N1 零模型（原始命中级）

Biomni 声明 B=1000 的 N1/N2 零模型，但 `results/layers/` 里没有交付任何比对表，
`data/` 也是空的 —— 零模型无法被任何第三方复算。本复现重跑 N1：

- 保组成洗牌 EBNA1（`np.random.default_rng(42+b)` 均匀置换），seeded 重扫，B=200
- 每次洗牌的命中数：**mean 48,096，sd 3,134，CV 0.065**

对照：
- 真实 seeded 59,108 / 洗牌 48,096 = **1.23×**
- 真实全量 139,363 / 洗牌 48,096 = **2.90×**

（Biomni 报告 N1 = 58,429、CV = 0.0008。我测到的是 48,096 与 CV 0.065 —— 命中数的
洗牌间波动比它报的大 80 倍，且中心值低 18%。它的 CV 数量级不可能来自"每次洗牌命中数"，
很可能量的是别的统计量。这是需要在正式稿里澄清的一处。）

---

## 3. 三个统计量并列，问题就暴露了

对每条真实命中，统计"零模型里有多少条比它更好的命中"（这是**计数**，不是 p 值）：

| 统计量 | 定义 | ==0 的命中 | < B 的命中 | 中位数 |
|---|---|---|---|---|
| `E_raw` | 洗牌命中 score ≥ 观测 score | 108 (0.08%) | 551 (0.40%) | 3,379,152 |
| `E_dens` | 洗牌 score/len ≥ 观测 | **0** | **0** | 6,666,925 |
| `E_band` | **同长度带内**洗牌 score ≥ 观测 | 420 (0.30%) | 953 (0.68%) | 1,477,461 |

读法：
- `E_raw = 0` 基本等价于"分数最高"，因为分数随比对长度单调增长 → 它奖励的是**最长**的比对，
  不是最有信息的比对。LORICRIN 291 aa / score 542 的 `E_raw=0`，但密度只有 1.86。
- `E_dens` 反过来被最短的比对淹没（6 aa / score 27 → 密度 4.5），**一条都过不了 null**。
- `E_band` 是唯一诚实的问法：「长度已经固定了，这个分数还异常吗？」

**关键：`E_band` 的 420 条 E=0 里，长度带内零模型完全无法产生同分数的比对 ——
但那是"洗牌把连续的 Gly-Ala 重复区打散"造成的结构性不可能，不是生物学信号。
所以对低复杂度驱动的命中，保组成洗牌**在原理上**就不能当零模型。**（这一点 Biomin 的
`docs/STATISTICS_INTERPRETATION.md` 只说了 E=0 是右删失，没说清它其实是结构性的 0。）

三个已发表基准落在哪：

| 基准 | score | 长度 | identity | E_band |
|---|---|---|---|---|
| MEF2B | 46 | 77 | 0.30 | 591 |
| ANO2 | 46 | 32 | 0.34 | 127,305 |
| HEPACAM | 43 | 33 | 0.45 | 158,321 |
| CRYAB | 32 | 32 | 0.38 | 271,264 |

即：**四个已知真阳性的比对层成绩，与保组成洗牌噪声不可区分。** 这与 Biomni 自己的
`STATISTICS_INTERPRETATION.md` §2（已验证真阳性 E_emp 671–11,146，不能用阈值截断）一致，
只是它的排序视图没有遵守这句话。

---

## 4. 去掉组成假象后还剩什么

用 Biomni 自己**定义过但从未用于排序**的过滤器，从 139,363 条往下筛：

| 过滤 | 剩余 | 占比 |
|---|---|---|
| 全部命中 | 139,363 | 100% |
| `aln_len ≥ 20` | 45,487 | 32.6% |
| 查询片段非 GAr 富集（gar_q < 0.50） | 26,453 | 19.0% |
| 靶标片段有序（entropy ≥ 2.8 且 maxfrac ≤ 0.30） | 20,604 | 14.8% |
| 靶标非 GA/SPR 重复（ga_t < 0.40 且 spr_t < 0.50） | 19,544 | 14.0% |
| **优于同长度带零模型（E_band < B）** | **8** | 0.006% |
| identity ≥ 0.35 | **1** | 0.001% |

**唯一幸存者：CBX8（Q9HC52），exhaustive，score 107，83 aa，identity 37.3%，
EBNA1 354–432 ↔ CBX8 162–244，gar_q 0.44，靶标 entropy 3.36，E_band 107。**

E_band=107/200 在 B=200 下只是勉强，B=1000 未必站得住。四个已发表基准**一个都没滤进来**。

---

## 5. 结论（当前证据支持的最强表述）

1. Biomni 的比对层可被独立复现，误差 0.04%。
2. 但它的排序主键 `E_emp(升序)` 在结构上等价于"最长比对优先"，六视图里 5 个的 top100
   有 90–100% 是 Gly-Ala 重复区 ↔ 胶原 Gly-X-Y 的组成假象。我复现出的 top-20 就是
   LORICRIN + 一堆 COL + ELN + HRNR + KRT，与它的视图 1 完全一致。
3. **在 EBNA1 这个特定查询上，序列比对层对分子模拟没有可检测的判别力**：
   已知真阳性与洗牌噪声同分布，去掉组成假象后只剩 1 条勉强边缘的 CBX8。
4. 这不是管线 bug，是方法边界：EBNA1 的序列被一段长 Gly-Ala 重复区主导，
   任何基于局部比对的打分都会被它带偏；而真正的 mimicry（若有）发生在结构/界面层。

因此 v1.3 报告里"潜在 mimicry 候选"的可信度，**完全来自注释层（文献 epitope 重叠、
PTM、HPA 表达），不来自比对层** —— 报告首页应明确这样写。

---

## 6. 复现命令

```bash
cd /workspace/PH-Mimic
gcc -O3 -march=native -pthread -o env/bin/mimic_scan src/mimic_scan.c

env/bin/mimic_scan --query data/cache/P03211.fasta --db data/cache/human_swissprot.fasta \
  --mode seeded     --out reports/seeded_raw.tsv
env/bin/mimic_scan --query data/cache/P03211.fasta --db data/cache/human_swissprot.fasta \
  --mode exhaustive --n-iter 4 --threads 8 --out reports/exhaustive_raw.tsv

env/venv/bin/python src/03_null_stats.py --real reports/seeded_raw.tsv,reports/exhaustive_raw.tsv \
  --out reports/null_stats_B200.tsv --B 200 --workers 10
env/venv/bin/python src/04_residual.py reports/null_stats_B200.tsv
```

## 7. 已知限制（写进正式稿之前必须处理）

- B=200 太粗，E=0 右删失在 1/B；要下 E<10 的判断需 B≥1000（约 2.5 min，可跑）。
- exhaustive 命中用的是**seeded** 零模型作参照分布，不是同模式零模型。
- `E_band` 的带边界（6/15/30/60/120/250/500）是任意的，需做敏感性分析。
- 洗牌零模型对连续低复杂度区**结构性无效**，需换成"局部保组成"零模型
  （窗口内洗牌 / 二肽保真 / 或 MMseqs2 的真实库 E-value）。
- identity ≥ 0.35 这个阈值没有文献依据，是本步为收口而设的探索性阈值。
