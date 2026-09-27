# 共线性分析与Ka/Ks分析操作文档

本文说明如何从基因组注释和蛋白质序列准备MCScanX输入、识别目标基因相关共线性关系、提取CDS并计算Ka/Ks。命令以Bash和`python3`为例；文件名、ID字段、参数及输出路径需按实际数据调整。开始前应固定参考基因组、注释、蛋白质和CDS文件的来源及版本，避免混用不同版本的数据。

## 一、分析设计与生物学前提

共线性分析用于识别不同基因组区域中基因顺序和邻接关系的保守性，可为基因组结构演化和基因复制历史提供证据。共线性伙伴不必然是直系同源基因；旁系同源关系、基因组复制、局部复制、基因丢失、注释差异以及组装质量均可能影响结果。种内分析和种间分析应分别明确比较对象、参考版本、ID体系及候选基因筛选规则；进行多物种或跨物种比较时，需确保不同物种的基因ID不会重名，并保留可追溯的物种前缀或映射表。

Ka/Ks（也常记作ω）表示非同义替换率Ka与同义替换率Ks之比。其解释依赖同源编码序列、密码子比对质量、替代模型、序列间分化程度和基因复制背景。Ka/Ks小于1通常与纯化选择相符，接近1通常与中性演化相符，大于1可提示正选择，但单个基因对的比值不能单独构成正选择证据。Ks较高时可能发生替换饱和；Ka或Ks接近0时，比值也可能不稳定。因此应结合原始Ka和Ks、比对质量、基因功能、物种关系及多基因层面的证据进行解释。

## 二、软件与目录准备

推荐在Linux、macOS或WSL中运行命令行步骤。Python脚本只使用标准库；流程需要Python 3、DIAMOND和TBtools-II。gffread用于在没有CDS FASTA时从基因组序列和注释提取CDS；MCScanX可由TBtools-II封装运行，也可单独使用命令行版本。安装后检查：

```bash
python3 --version
diamond version
gffread --version    # 仅在需要时检查
```

建议将文件按以下方式组织，并将真实输入、临时文件和结果保存在仓库中相应目录：

```text
data/raw/genome.gff
data/raw/protein.fasta
data/raw/all.cds.fa       # 已有全量CDS时使用
data/raw/genome.fasta     # 仅在需要用gffread生成CDS时使用
data/intermediate/
results/mcscanx/
results/kaks/
```

## 三、输入文件与ID对应关系

GFF/GFF3必须能从选定的feature及属性字段中提取与蛋白质FASTA相对应的ID。蛋白质FASTA标题默认以`>`之后、空格之前的第一个字段作为蛋白质ID；CDS FASTA的标题也必须能与候选基因对中的ID精确对应。若使用基因ID、转录本ID和蛋白质ID中的不同层级，须先建立明确映射，不能仅凭相似字符串猜测对应关系。一个基因有多个转录本时，应预先规定代表转录本选择策略，避免在共线性和Ka/Ks步骤混用不同转录本。

输入要求如下：

| 文件 | 要求 |
|---|---|
| GFF/GFF3 | feature类型及属性字段须能稳定提取目标ID和坐标；坐标须与所用基因组版本一致 |
| 蛋白质FASTA | 每条序列有唯一、可追溯的ID，且与注释中提取的ID相符 |
| CDS FASTA | 用于Ka/Ks；CDS ID须与候选基因对中的ID完全一致 |
| 基因组FASTA | 仅在需要用gffread生成CDS时必需，须与GFF/GFF3版本相匹配 |

先抽查GFF中的CDS记录和蛋白质FASTA标题：

```bash
awk -F '\t' '$3 == "CDS" {print; n++; if (n == 3) exit}' data/raw/genome.gff
grep -m 3 '^>' data/raw/protein.fasta
```

脚本默认从GFF/GFF3的`CDS`行读取`Protein_Accession`属性。若实际注释使用`protein_id`，或蛋白质ID对应于`mRNA`行的`ID`，应按真实文件指定参数。不要直接照搬示例；还须确认所选feature能产生适用于MCScanX的基因位置，并检查生成结果中每个ID是否唯一、坐标是否合理。若CDS由多个外显子记录组成，应确认脚本按预期汇总同一基因或蛋白质的坐标，而不是将每个CDS片段误当成独立基因。

## 四、生成MCScanX输入

若GFF使用脚本默认的`CDS`和`Protein_Accession`：

```bash
python3 scripts/prepare_mcscanx_inputs.py \
  --gff data/raw/genome.gff \
  --protein data/raw/protein.fasta \
  --prefix TK \
  --outdir data/intermediate
```

若蛋白质ID位于GFF的`protein_id`属性：

```bash
python3 scripts/prepare_mcscanx_inputs.py \
  --gff data/raw/genome.gff \
  --protein data/raw/protein.fasta \
  --prefix TK \
  --outdir data/intermediate \
  --feature_type CDS \
  --protein_attr protein_id
```

若ID对应于`mRNA`行的`ID`属性，则应按注释格式调整为`--feature_type mRNA --protein_attr ID`。生成文件通常为`data/intermediate/TK.clean.protein.fasta`和`data/intermediate/TK.gff`。后者为MCScanX四列格式：染色体或scaffold、基因或蛋白质ID、起始坐标、终止坐标。

脚本会打印蛋白质FASTA记录数、GFF中匹配到坐标的ID数及写入行数。写入行数为0，或明显低于预期有效编码基因数时，应先检查feature类型、属性字段、ID版本号、转录本后缀和FASTA标题中的描述字段，不要继续运行DIAMOND。还应检查`TK.gff`中ID唯一性、染色体命名和坐标范围。种间分析中需确保各物种ID唯一，且染色体或scaffold命名可区分物种。

## 五、DIAMOND蛋白质全对全比对

在中间文件目录创建数据库并运行`blastp`：

```bash
cd data/intermediate

diamond makedb \
  --in TK.clean.protein.fasta \
  --db TK

diamond blastp \
  --db TK.dmnd \
  --query TK.clean.protein.fasta \
  --out TK.blast \
  --evalue 1e-5 \
  --max-target-seqs 5 \
  --outfmt 6

cd ../..
```

`--evalue 1e-5`和每条查询最多保留5个匹配仅为示例参数，不是普适标准。对于重复基因丰富的基因组或需要保留更多候选同源关系的研究，过低的`--max-target-seqs`可能遗漏可用于后续分析的匹配，应根据研究目标调整，并在方法中记录E-value、匹配数量上限、输出字段和软件版本。生成的主要文件为`data/intermediate/TK.dmnd`和`data/intermediate/TK.blast`。

## 六、运行MCScanX并检查共线性结果

在TBtools-II中打开`Quick Run MCScanX Wrapper`，设置比对文件为`data/intermediate/TK.blast`、GFF文件为`data/intermediate/TK.gff`，输出目录为`results/mcscanx`。核心输出应包括`results/mcscanx/TK.collinearity`。若直接运行命令行版MCScanX，应遵循该版本要求的文件前缀、目录和参数设置，并记录软件版本及命令。

运行后应确认`.collinearity`文件非空，基因ID与`TK.gff`及DIAMOND比对文件一致，并抽查若干共线性区块及其基因对。若ID不一致、输出异常少或无有效区块，应回查蛋白质FASTA、ID映射、GFF四列坐标及DIAMOND筛选参数，而不是直接解释生物学结果。

## 七、筛选目标基因相关的共线性基因对

复制目标基因示例文件并编辑，每行填写一个目标基因或蛋白质ID：

```bash
cp configs/targets.example.txt configs/targets.txt
```

```text
GENE000001
GENE000002
GENE000003
```

运行筛选脚本：

```bash
python3 scripts/filter_collinearity_targets.py \
  --collinearity results/mcscanx/TK.collinearity \
  --gff data/intermediate/TK.gff \
  --targets configs/targets.txt \
  --out_prefix results/mcscanx/target_collinearity
```

脚本生成`target_collinearity_pairs.tsv`和`target_collinearity_summary.tsv`。前者记录目标基因、共线性伙伴、区块和坐标；后者汇总目标基因是否命中、伙伴基因数、区块数及同染色体或异染色体关系。`found_in_collinearity = NO`时，先确认目标列表、`TK.gff`和`.collinearity`使用的是同一ID层级和命名规则。该脚本按MCScanX标准`.collinearity`文本格式解析；若软件版本或输出格式不同，应先抽查解析结果。

## 八、准备Ka/Ks基因对列表

复制示例并编辑：

```bash
cp configs/kaks_pairs.example.tsv configs/kaks_pairs.tsv
```

`kaks_pairs.tsv`每行填写一对基因ID，以Tab分隔且不添加表头。示例中的`<Tab>`仅表示制表符，实际文件中不要输入尖括号及单词`Tab`：

```text
GENE000001<Tab>GENE000002
GENE000003<Tab>GENE000004
```

生成去重后的ID列表：

```bash
awk '{print $1; print $2}' configs/kaks_pairs.tsv | sort -u > configs/kaks_ids.txt
```

确认每一对序列的来源和关系，并记录基因对筛选标准。共线性只能作为关系判断的上下文证据之一；种内复制基因对与跨物种基因对的演化解释不同，不应混为同一类比较。

## 九、提取候选CDS

若已有全量CDS FASTA：

```bash
mkdir -p results/kaks
python3 scripts/extract_cds_by_ids.py \
  --cds data/raw/all.cds.fa \
  --ids configs/kaks_ids.txt \
  --out results/kaks/target.cds.fa
```

若没有CDS FASTA，但有与注释版本一致的基因组FASTA和GFF/GFF3，可先用gffread提取：

```bash
gffread data/raw/genome.gff \
  -g data/raw/genome.fasta \
  -x data/intermediate/all.cds.fa

python3 scripts/extract_cds_by_ids.py \
  --cds data/intermediate/all.cds.fa \
  --ids configs/kaks_ids.txt \
  --out results/kaks/target.cds.fa
```

gffread输出的ID可能是转录本ID，而候选基因对使用的可能是基因ID或蛋白质ID；运行前抽查FASTA标题及ID映射：

```bash
grep -m 5 '^>' data/intermediate/all.cds.fa
head configs/kaks_ids.txt
```

提取脚本会报告找到和缺失的ID。`Found`应与请求ID数量相符；存在`Missing`时，应回查ID映射、版本号和代表转录本选择规则，不要仅通过模糊字符串替换强行对应。

## 十、CDS质量检查

```bash
python3 scripts/validate_cds_lengths.py \
  --cds results/kaks/target.cds.fa \
  > results/kaks/cds_validation.tsv
```

检查结果包括CDS长度、长度是否为3的倍数、是否以ATG起始、是否以TAA/TAG/TGA终止及是否含ATCGN以外字符。这些项目属于基础质量控制提示，并非自动删除规则：部分合法基因模型可能不含标准起始或终止密码子；模糊碱基、遗传密码表、叶绿体或线粒体编码规则也可能影响判断。Ka/Ks分析前还须确认两条序列对应正确的编码转录本、阅读框一致，且密码子比对质量可接受。

## 十一、计算Ka/Ks

在TBtools-II中打开`Simple Ka/Ks Calculator (NG)`，设置CDS文件为`results/kaks/target.cds.fa`、基因对文件为`configs/kaks_pairs.tsv`，输出文件为`results/kaks/kaks.result.tsv`。完成后检查Gene 1、Gene 2、Ka、Ks、Ka/Ks以及失败、NaN或无穷值，并记录TBtools-II版本、方法和参数。确保候选序列按所选方法进行合理的密码子比对；未经检查的CDS逐碱基比较可能因插缺、移码或错误转录本导致无效结果。

## 十二、结果判读与报告

| Ka/Ks或ω | 常见初步解释 | 注意事项 |
|---|---|---|
| `< 1` | 与纯化选择相符 | 需排除序列错配、错误CDS及Ks估算问题 |
| `≈ 1` | 与中性演化相符 | 仅是简化解释，受模型和统计不确定性影响 |
| `> 1` | 可能提示正选择或功能分化 | 单个基因对不足以证明适应性进化，应结合序列、模型和多基因证据 |
| `NaN`或`Inf` | 当前数据或模型下无法稳定估算 | 不代表纯化选择，也不应直接替换为0 |

Ka/Ks是针对一对同源编码序列的估计值，不是单个基因固有属性。Ks过高可能发生替换饱和；Ka或Ks接近0时比值可能不稳定。方法报告应至少说明参考基因组、注释与蛋白质序列版本、ID映射规则、DIAMOND/MCScanX/TBtools-II/gffread版本、DIAMOND参数、输入和匹配数量、共线性筛选规则、代表转录本策略、CDS质量控制方法、Ka/Ks计算方法，以及NaN、Inf和饱和Ks的处理原则。

## 十三、常见问题

**`Rows written`为0。** 通常表示`--feature_type`或`--protein_attr`与GFF实际内容不一致。检查GFF第三列feature类型和第九列属性字段，并确认提取的ID与蛋白质FASTA标题一致。

**只有少量蛋白质写入`TK.gff`。** 检查蛋白质FASTA与注释中的ID是否只部分匹配，重点核对版本号、前缀、转录本后缀、空格后的描述及ID层级。

**TBtools-II或MCScanX提示ID不一致。** 确认`TK.clean.protein.fasta`标题、`TK.blast`前两列和`TK.gff`第二列采用兼容且唯一的ID。

**目标基因均未命中共线性结果。** 确认目标列表使用的ID层级与`TK.gff`及`.collinearity`完全一致，并确认目标ID存在于输入数据中。

**CDS提取出现Missing。** 回到注释和FASTA标题建立明确的gene、transcript、protein映射，并检查是否使用了正确的代表转录本。

**Ka/Ks结果为NaN或Inf。** 可能与有效替换位点不足、序列相同或近似相同、Ka或Ks为0、阅读框错误、密码子比对异常或模型不适用有关。应检查原始序列、比对和模型，不要将NaN解释为某类选择压力或直接改写为0。

## 十四、数据安全与提交前检查

`.gitignore`默认排除FASTA、GFF/GTF、GenBank文件、DIAMOND数据库与比对结果、MCScanX结果、常见表格、图片、日志和报告，以及`data/raw/`、`data/intermediate/`和`results/`中的实际内容。忽略规则不能替代提交前审查；上传前运行：

```bash
git status
```

确认没有真实序列、未公开数据、样本信息或不宜公开的分析结果被加入提交。
