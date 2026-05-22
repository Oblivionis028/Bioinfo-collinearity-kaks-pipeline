<h1 align="center">🧬 共线性分析与 Ka/Ks 分析流程模板</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.x-3776AB?style=flat-square&logo=python&logoColor=white" />
  <img src="https://img.shields.io/badge/Shell-Bash-4EAA25?style=flat-square&logo=gnubash&logoColor=white" />
  <img src="https://img.shields.io/badge/MCScanX-Collinearity-2E8B57?style=flat-square" />
  <img src="https://img.shields.io/badge/Ka%2FKs-Analysis-DC143C?style=flat-square" />
  <img src="https://img.shields.io/badge/Data-Not%20Included-lightgrey?style=flat-square" />
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" />
</p>

一个可重复使用的 **种内/种间共线性分析** 与 **Ka/Ks 分析**流程模板。

本仓库用于整理从 `GFF/GFF3` 注释文件和 `protein FASTA` 蛋白序列文件出发，完成 **MCScanX 共线性分析**、**目标基因筛选** 以及 **Ka/Ks 选择压力分析** 的标准流程。

> 本仓库仅保存分析流程、辅助脚本和示例配置，不包含任何真实基因组数据、注释文件、序列文件或分析结果。

## ✨ 功能特点

- 🧬 从 GFF/GFF3 和蛋白 FASTA 生成 MCScanX 输入文件
- ⚡ 使用 DIAMOND 进行全基因组蛋白自比对
- 🔗 使用 MCScanX / TBtools 进行共线性分析
- 🎯 从共线性结果中筛选目标基因相关的共线性基因对
- 🧪 提取候选基因对的 CDS 序列用于 Ka/Ks 分析
- 📊 根据 Ka/Ks 结果判断基因对受到的选择压力

## 仓库结构

```text
.
├── scripts/
│   ├── prepare_mcscanx_inputs.py
│   ├── filter_collinearity_targets.py
│   ├── extract_cds_by_ids.py
│   └── validate_cds_lengths.py
├── configs/
│   ├── targets.example.txt
│   └── kaks_pairs.example.tsv
├── data/
│   ├── raw/              # 本地放原始数据；不会上传 GitHub
│   └── intermediate/     # 中间文件；不会上传 GitHub
├── results/
│   ├── mcscanx/          # MCScanX 结果；不会上传 GitHub
│   └── kaks/             # Ka/Ks 结果；不会上传 GitHub
├── docs/
│   └── workflow.zh-CN.md
├── .gitignore
├── LICENSE
└── README.md
```

## 数据隐私

本仓库已经配置 `.gitignore`，默认不会上传以下敏感或大文件：

- `*.fa`
- `*.fasta`
- `*.fna`
- `*.gff`
- `*.gff3`
- `*.gtf`
- `*.blast`
- `*.dmnd`
- `*.collinearity`
- `data/raw/*`
- `data/intermediate/*`
- `results/*`

上传 GitHub 前仍建议运行：

```bash
git status
```

确认没有真实序列、真实注释或真实分析结果被加入 Git。

## 快速流程

### 1. 准备原始文件

把原始文件放在本地，不提交到 GitHub：

```text
data/raw/genome.gff
data/raw/protein.fasta
data/raw/genome.fasta   # 仅 Ka/Ks 需要从 genome + gff 提取 CDS 时使用
```

### 2. 生成 MCScanX 输入

```bash
python3 scripts/prepare_mcscanx_inputs.py \
  --gff data/raw/genome.gff \
  --protein data/raw/protein.fasta \
  --prefix TK \
  --outdir data/intermediate
```

输出：

```text
data/intermediate/TK.clean.protein.fasta
data/intermediate/TK.gff
```

其中 `TK.gff` 是 MCScanX 4 列格式：

```text
chromosome/scaffold    gene_id/protein_id    start    end
```

### 3. DIAMOND 蛋白自比对

```bash
cd data/intermediate

diamond makedb --in TK.clean.protein.fasta -d TK
diamond blastp -d TK.dmnd -q TK.clean.protein.fasta -o TK.blast -e 1e-5 -k 5 -f 6
```

输出：

```text
TK.blast
```

### 4. MCScanX / TBtools 共线性分析

打开 TBtools：

```text
Quick Run MCScanX Wrapper
```

填入：

```text
Input .blast File: data/intermediate/TK.blast
Input .gff File:   data/intermediate/TK.gff
Output Directory:  results/mcscanx
```

核心输出：

```text
results/mcscanx/TK.collinearity
```

### 5. 筛选目标基因相关共线性基因对

准备目标基因列表：

```text
configs/targets.txt
```

一行一个目标基因/蛋白 ID。

运行：

```bash
python3 scripts/filter_collinearity_targets.py \
  --collinearity results/mcscanx/TK.collinearity \
  --gff data/intermediate/TK.gff \
  --targets configs/targets.txt \
  --out_prefix results/mcscanx/target_collinearity
```

输出：

```text
results/mcscanx/target_collinearity_pairs.tsv
results/mcscanx/target_collinearity_summary.tsv
```

### 6. Ka/Ks 候选基因对

准备候选基因对：

```text
configs/kaks_pairs.tsv
```

格式：

```text
geneA    geneB
```

中间是 Tab。

### 7. 提取 CDS

如果已经有全基因组 CDS FASTA：

```bash
python3 scripts/extract_cds_by_ids.py \
  --cds data/raw/all.cds.fa \
  --ids configs/kaks_ids.txt \
  --out results/kaks/target.cds.fa
```

如果没有 CDS FASTA，但有 genome FASTA + GFF，可先用 `gffread` 提取：

```bash
gffread data/raw/genome.gff -g data/raw/genome.fasta -x data/intermediate/all.cds.fa
```

### 8. 检查 CDS 长度

```bash
python3 scripts/validate_cds_lengths.py --cds results/kaks/target.cds.fa
```

### 9. 计算 Ka/Ks

打开 TBtools：

```text
Simple Ka/Ks Calculator (NG)
```

填入：

```text
Input CDS File:      results/kaks/target.cds.fa
Input GenePair File: configs/kaks_pairs.tsv
Output Table File:   results/kaks/kaks.result.tsv
```

## 结果解释

常见解释：

| Ka/Ks | 解释 |
|---|---|
| `< 1` | 纯化选择，功能较保守 |
| `≈ 1` | 中性进化 |
| `> 1` | 正选择，可能存在功能分化 |
| `NaN` | Ka 和 Ks 不能有效估计，通常不能用于选择压力判断 |

## 依赖

- Python 3
- DIAMOND
- TBtools-II
- 可选：gffread

## 注意

本模板不包含任何真实研究数据。请勿提交真实序列、真实注释、真实 `.blast`、真实 `.collinearity` 或真实 Ka/Ks 结果。
