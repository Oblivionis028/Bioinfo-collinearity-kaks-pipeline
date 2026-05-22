# Workflow：共线性分析与 Ka/Ks 分析

## 1. 输入文件说明

最少需要：

```text
protein.fasta
annotation.gff / annotation.gff3
```

用于共线性分析。

Ka/Ks 还需要：

```text
cds.fasta
```

如果没有 CDS FASTA，但有：

```text
genome.fasta
annotation.gff
```

可以用 `gffread` 提取 CDS。

## 2. 共线性分析核心逻辑

```text
protein.fasta
↓
DIAMOND / BLASTP all-vs-all
↓
blast tabular result
↓
MCScanX + gene position file
↓
collinearity blocks
```

## 3. Ka/Ks 分析核心逻辑

```text
候选同源基因对
↓
提取双方 CDS
↓
进行密码子层面的替换率计算
↓
获得 Ka、Ks、Ka/Ks
```

Ka/Ks 不是对单个基因算，而是对一对同源基因算。

## 4. 推荐结果记录

建议记录：

- 全基因组输入基因数
- 共线性基因数
- 共线性基因比例
- 目标基因命中数量
- 去重后的目标相关共线性基因对
- Ka、Ks、Ka/Ks
- 对 NaN 结果的说明

## 5. 常见问题

### MCScanX 只识别出很少基因

检查 4 列 GFF 的列顺序是否正确：

```text
chromosome/scaffold    gene_id    start    end
```

不要写成：

```text
gene_id    chromosome/scaffold    start    end
```

### TBtools One Step MCScanX 报 ID 不一致

建议改用标准流程：

```text
DIAMOND 生成 blast
+
4列 GFF
+
Quick Run MCScanX Wrapper
```

### Ka/Ks 出现 NaN

如果 Ka = 0 且 Ks = 0，则 Ka/Ks 为 0/0，无法估计。该基因对不适合用于选择压力判断。
