# 共线性分析与Ka/Ks分析操作文档

本文按照“准备文件→检查ID→生成MCScanX输入→蛋白质序列比对→识别共线性区块→筛选目标基因对→提取CDS→质量检查→计算Ka/Ks→解释结果”的顺序，说明如何完成一轮可追溯的分析。文档分为两条环境路线：第2至第14节为Ubuntu/Linux终端（Bash）说明；第15节提供Windows原生环境（PowerShell）的对应步骤。WSL实际运行的是Linux环境，应遵循Linux路线，不能把Bash命令直接粘贴到PowerShell。TBtools-II图形界面步骤可用于两条路线，但选择文件时应使用当前操作系统的路径格式。示例文件名和参数用于演示，必须按实际物种、注释格式和研究问题调整。本文中的脚本路径和参数对应本仓库当前版本；若本地脚本版本不同，应先查看对应脚本的帮助信息。

共线性分析的最低输入是蛋白质FASTA和GFF/GFF3注释；Ka/Ks分析还需要CDS FASTA。如果没有CDS FASTA，但有与注释版本完全匹配的基因组FASTA及GFF/GFF3，可尝试用gffread提取。完整路线为：蛋白质序列→DIAMOND蛋白质全对全比对→比对表格和基因位置文件→MCScanX识别共线性区块→筛选目标基因相关候选基因对→提取双方CDS→核对CDS和密码子比对质量→计算Ka、Ks及Ka/Ks。仓库脚本负责输入准备、目标筛选、CDS提取和基础检查；DIAMOND、MCScanX及Ka/Ks计算仍由外部软件完成，因此不能仅凭“脚本运行成功”判断生物学结果正确。

开始前先完成四项核对：第一，确认蛋白质、CDS、GFF/GFF3和基因组FASTA（如使用）来自同一物种和同一基因组注释版本；第二，确定本次是种内还是种间分析，以及分析对象和目标基因；第三，确认GFF/GFF3中哪一列属性对应蛋白质FASTA和CDS FASTA的ID；第四，确定一个基因有多个转录本时采用何种代表转录本规则，并在后续步骤中始终遵循该规则。

## 一、分析设计与生物学前提

共线性分析用于识别不同基因组区域中基因顺序和邻接关系的保守性，可为基因组结构演化和基因复制历史提供证据。共线性伙伴不必然是直系同源基因；旁系同源关系、基因组复制、局部复制、基因丢失、注释差异以及组装质量均可能影响结果。种内分析和种间分析应分别明确比较对象、参考版本、ID体系及候选基因筛选规则；进行多物种或跨物种比较时，需确保不同物种的基因ID不会重名，并保留可追溯的物种前缀或映射表。

本操作文档的单物种命令示例用于种内共线性分析。种间分析需要A、B两个物种各自的蛋白质FASTA和位置文件；两套位置文件合并后，第一列中的染色体/scaffold名称及第二列基因ID都必须能区分物种。若A、B都使用`Chr1`或都存在同名基因ID，需在分析前建立带物种前缀的唯一名称，并同步更新所有相关文件及映射记录。准备脚本的`--prefix`只更改输出文件名，不会自动修改序列ID或染色体名。

种间分析也不建议不加说明地将两物种蛋白质拼成一个FASTA、对合并集合只做一次全对全搜索。MCScanX原始说明允许合并多物种的种内与种间比对表；其论文进一步建议，多物种分析可将种内与物种间分别产生的BLASTP结果合并，并对保留匹配设置best-hit筛选，以避免近缘物种匹配过多造成偏倚。下面第5节提供A↔A、B↔B、A↔B、B↔A四个DIAMOND搜索的示例。MCScanX的不同入口对位置文件格式可能有差异，因此此处仍按本仓库脚本与TBtools-II`Quick Run MCScanX Wrapper`使用四列位置文件；切换到其他MCScanX入口前须复核格式。

分别运行A、B的输入准备步骤后，用以下命令合并两份无表头四列位置文件。`cat`只适用于两份文件的列顺序完全一致，且ID和染色体/scaffold名称都已全局唯一的情形：

```bash
cat data/intermediate/A.gff data/intermediate/B.gff \
  > data/intermediate/AB.gff
```

跨物种Ka/Ks筛选出的候选对也必须使用与两物种蛋白质/CDS完全一致的带前缀ID；不能只给`AB.gff`中的ID加前缀。关于MCScanX多物种比对输入的说明见[MCScanX原始README](https://github.com/wyp1125/MCScanX#mcscanx)及[MCScanX论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336/)。

Ka/Ks（也常记作ω）表示非同义替换率Ka与同义替换率Ks之比。其解释依赖同源编码序列、密码子比对质量、替代模型、序列间分化程度和基因复制背景。Ka/Ks小于1通常与纯化选择相符，接近1通常与中性演化相符，大于1可提示正选择，但单个基因对的比值不能单独构成正选择证据。Ks较高时可能发生替换饱和；Ka或Ks接近0时，比值也可能不稳定。因此应结合原始Ka和Ks、比对质量、基因功能、物种关系及多基因层面的证据进行解释。

## 二、软件与目录准备

本节给出Linux环境的安装和目录准备。Linux/WSL适合完整命令行流程；Windows原生安装方法及其兼容性边界见第15节。以下安装示例使用Conda环境隔离依赖，不会改动系统Python。若尚未安装Conda，可按[Miniforge官方说明](https://github.com/conda-forge/miniforge)安装；不要将不同教程里的多个环境管理器混装到同一个环境中。本仓库的Python脚本仅使用Python标准库，不需要通过`pip`额外安装包。DIAMOND负责蛋白质比对，MCScanX负责根据比对和基因位置识别共线性区块；TBtools-II可提供图形界面封装及Ka/Ks计算；只有需要从基因组FASTA提取CDS时才需要gffread。

如果在Ubuntu/Linux x86_64上尚未安装Conda，可先在浏览器打开Miniforge官方仓库的Releases页面，下载最新的`Miniforge3-Linux-x86_64.sh`安装文件到“下载”目录。在Linux终端中运行：

```bash
cd ~/Downloads
bash Miniforge3-Linux-x86_64.sh
```

安装程序会逐步询问许可和安装位置：阅读许可后输入`yes`接受；安装位置可使用默认的`~/miniforge3`；当程序询问是否初始化Conda时选择`yes`。安装结束后关闭并重新打开终端，再输入`conda --version`确认命令可用。若电脑不是x86_64架构，不要下载上述文件名，应在官方页面选择与`uname -m`结果相符的安装包。若`conda activate`提示尚未初始化，可在Bash终端执行`conda init bash`，再关闭并重新打开终端。本文第2至第14节的命令均为Bash语法；Windows用户请切换至第15节的PowerShell步骤。

在终端中创建专用环境并安装命令行工具：

```bash
conda create -n collinearity python=3.11 diamond gffread mcscanx \
  --channel conda-forge --channel bioconda --strict-channel-priority
conda activate collinearity
```

若当前只做已有CDS FASTA的分析，可不安装gffread；若通过TBtools-II图形界面运行MCScanX，仍建议安装MCScanX命令行程序以便必要时独立检查版本或排错。TBtools-II应从其[官方仓库及发布页面](https://github.com/CJ-Chen/TBtools-II)单独安装：Windows使用与系统匹配的安装包；Linux/macOS使用官方发布的跨平台包，按发布说明解压并启动；第一次启动后确认主窗口能正常打开，并在软件界面中记录版本号。若跨平台包无法启动，再根据该发布版本的要求检查Java环境，不要为解决启动问题而随意安装来源不明的旧版Java。Conda仓库和软件版本会更新，因此正式分析时应保存版本信息，不要只记录“使用最新版”。

激活环境后检查：

```bash
python3 --version
diamond version
gffread --version    # 仅在使用gffread时检查
command -v MCScanX   # 仅在使用命令行MCScanX时检查
conda list
```

记录上述输出中的Python、DIAMOND、gffread、MCScanX和TBtools-II版本。若希望保留当前Linux平台的精确Conda软件构建清单，可在仓库根目录运行`conda list --explicit > conda-linux-64.explicit.txt`；该清单通常与操作系统及CPU平台有关，不能直接当作跨平台安装文件。另可用`conda env export --from-history > environment.yml`记录主动安装的软件需求。保存清单时应确认其位于预期目录，并在提交仓库前检查是否包含不应公开的路径或信息。

如尚未取得仓库，可在终端中克隆并进入项目目录；已有本地仓库则直接进入该目录：

```bash
git clone https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline.git
cd Bioinfo-collinearity-kaks-pipeline
```

将目录结构准备为如下形式。真实输入保存在`data/raw/`，自动生成或临时文件保存在`data/intermediate/`，正式分析结果保存在`results/`；目标列表和候选基因对也建议放在`data/intermediate/`，避免把真实研究配置覆盖仓库内的示例文件。

```text
data/raw/genome.gff
data/raw/protein.fasta
data/raw/all.cds.fa       # 已有全量CDS时使用
data/raw/genome.fasta     # 仅在需要用gffread生成CDS时使用
data/intermediate/
results/mcscanx/
results/kaks/
```

创建本次运行所需的目录：

```bash
mkdir -p data/raw data/intermediate results/mcscanx results/kaks
```

后续命令均假设当前工作目录为仓库根目录，可用`pwd`确认。若终端提示找不到`data/raw/`、`scripts/`等路径，通常是尚未进入仓库根目录，不要通过随意改写路径掩盖这个问题。

## 三、输入文件与ID对应关系

分析中最容易导致“软件成功运行但结果极少或为空”的问题，是不同文件使用了不同层级的ID。基因ID代表基因座，转录本ID代表某个剪接异构体，蛋白质ID代表翻译产物；三者常有不同命名。例如，蛋白质FASTA可能使用`XP_...`，GFF属性中同时出现基因ID、转录本ID和蛋白质ID，而gffread生成的CDS FASTA可能以转录本ID作为标题。必须根据注释发布方提供的字段定义建立明确映射，不能只因ID前后缀相似就自行删改。

本仓库脚本读取FASTA标题中`>`之后、空格之前的第一个字段作为序列ID。例如标题`>PROT001 hypothetical protein`的ID是`PROT001`。ID应唯一且在相关文件中精确一致。若FASTA标题中含有版本号，GFF注释或配对文件也必须采用同一形式。多物种分析时，建议在物种的基因、转录本及蛋白质ID前统一加入不含空格的物种前缀，例如`SpeciesA_gene001`，并确保蛋白质序列、位置文件、CDS序列及基因对列表采用一致ID；必要时也为染色体或scaffold名加物种前缀。不要只修改其中一个文件。

一个基因有多个转录本时，必须在全流程中采用一致策略，例如使用数据库指定的代表转录本，或按研究方案筛选具有完整CDS的转录本。最长转录本并不自动等同于最适合分析的代表转录本。对Ka/Ks而言，配对双方应是经研究设计筛选的同源编码序列；错误拼接不同异构体会改变阅读框或替换率估计。

输入要求如下：

| 文件 | 要求 |
|---|---|
| GFF/GFF3 | feature类型及属性字段须能稳定提取目标ID和坐标；坐标须与所用基因组版本一致 |
| 蛋白质FASTA | 每条序列有唯一、可追溯的ID，且与注释中提取的ID相符 |
| CDS FASTA | 用于Ka/Ks；CDS ID须与候选基因对中的ID完全一致 |
| 基因组FASTA | 仅在需要用gffread生成CDS时必需，须与GFF/GFF3版本相匹配 |

先确认文件确实存在且非空，并查看FASTA标题及注释内容。以下命令从仓库根目录执行；若数据文件名不同，替换成实际路径：

```bash
ls -lh data/raw/protein.fasta data/raw/genome.gff
grep -m 5 '^>' data/raw/protein.fasta
awk -F '\t' '$3 == "CDS" {print $1, $3, $4, $5, $9; n++; if (n == 5) exit}' data/raw/genome.gff
```

GFF/GFF3通常为Tab分隔的九列文件，第三列是feature类型，第四和第五列是坐标，第九列是属性。仓库的`prepare_mcscanx_inputs.py`默认只读取第三列恰好为`CDS`的记录，并从第九列读取`Protein_Accession=...`形式的键值属性；它能按`key=value`形式查找指定属性，不应假定它能直接解析GTF常用的`key "value"`格式。若实际字段是`protein_id=...`，后续命令需指定`--protein_attr protein_id`；若序列ID对应`mRNA`行中的`ID=...`，则需结合真实注释结构确认应选用何种feature和属性。若GFF压缩为`.gz`，应先解压或确认脚本支持后再处理。

实际查看时不要只看GFF的第一行，因为开头可能是注释行或其他feature。应找到一条与目标蛋白相对应的CDS记录，逐列确认：第三列是`CDS`，第九列中确实包含选定的属性键，等号后面的值与FASTA标题首字段相同。例如，若GFF属性为`ID=cds-XP_0001;Parent=rna-XP_0001;protein_id=XP_0001`，而蛋白质FASTA标题是`>XP_0001 predicted protein`，则应尝试`--protein_attr protein_id`；如果选定的属性值是`rna-XP_0001`，而FASTA ID是`XP_0001`，二者并不匹配，需查明注释说明或建立映射。可用以下命令只抽查含指定属性的CDS行：

```bash
awk -F '\t' '$3 == "CDS" && /Protein_Accession=/ {print $1, $3, $4, $5, $9; n++; if (n == 5) exit}' data/raw/genome.gff
awk -F '\t' '$3 == "CDS" && /protein_id=/ {print $1, $3, $4, $5, $9; n++; if (n == 5) exit}' data/raw/genome.gff
grep -m 5 '^>' data/raw/protein.fasta
```

如果第二条命令没有输出，说明GFF的CDS行中可能没有`protein_id=`这个属性，不能直接使用该参数；应根据注释文件中真实存在的属性键替换。抽查多个记录以确认这不是个别例外。

GFF/GFF3中的`seqid`（第一列）应与基因组FASTA的序列标题对应；否则gffread可能无法定位外显子或CDS片段。抽查两者命名：

```bash
grep -m 5 '^>' data/raw/genome.fasta    # 仅在准备使用gffread时检查
awk -F '\t' 'BEGIN {OFS="\t"} !/^#/ && NF >= 9 {print $1; n++; if (n == 5) exit}' data/raw/genome.gff
```

以上抽查若显示注释ID为`gene001.1`、蛋白质FASTA ID为`gene001.1.t1`、CDS标题为`transcript001`等不同格式，应在继续之前查阅注释文件说明或建立正式映射表。不要靠删除末尾字符、正则替换或模糊匹配批量“修复”ID，因为这可能把不同基因合并成同一个ID。

## 四、生成MCScanX输入

本步骤将蛋白质FASTA与GFF/GFF3坐标按ID匹配，生成清理后的蛋白质FASTA和MCScanX位置文件。先根据上一节的抽查结果选定正确的feature类型和ID属性，不确定时不要直接复制默认参数。若GFF第三列为`CDS`且第九列属性中存在`Protein_Accession=...`，运行：

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

若ID对应于`mRNA`行的`ID`属性，则可将参数调整为`--feature_type mRNA --protein_attr ID`，但前提是该ID确实与蛋白质FASTA标题一致，并且该feature的坐标可代表目标转录本。`--prefix TK`中的`TK`仅用于生成文件的命名前缀，不会自动给FASTA序列ID添加物种前缀；种间分析的ID去重必须提前在各输入文件中统一完成。

若做两物种比较，分别为A和B各运行一次准备脚本，并让输出文件名前缀不同。以下示例假设原始文件分别放在`data/raw/A.gff`、`data/raw/A.protein.fasta`、`data/raw/B.gff`和`data/raw/B.protein.fasta`，且两套文件内部ID已能对应：

```bash
python3 scripts/prepare_mcscanx_inputs.py \
  --gff data/raw/A.gff --protein data/raw/A.protein.fasta \
  --prefix A --outdir data/intermediate
python3 scripts/prepare_mcscanx_inputs.py \
  --gff data/raw/B.gff --protein data/raw/B.protein.fasta \
  --prefix B --outdir data/intermediate
```

两次运行后应看到`A.clean.protein.fasta`、`A.gff`、`B.clean.protein.fasta`和`B.gff`。`A`、`B`只是输出文件前缀，不会更改FASTA ID；如果两物种的ID重名，必须先在输入数据层面统一重命名，并确保重命名规则同步应用于对应GFF属性、蛋白质FASTA和CDS FASTA。若A、B的染色体名称重名，也需在合并位置文件前区分第一列的名称；如果之后通过gffread提取CDS，基因组FASTA中的序列标题亦须与各自GFF第一列保持对应。

生成文件通常为`data/intermediate/TK.clean.protein.fasta`和`data/intermediate/TK.gff`。本流程约定`TK.gff`为Tab分隔的四列位置文件，列序是染色体或scaffold、基因或蛋白质ID、起始坐标、终止坐标；以下仅为列顺序示意，不能把表头写入实际文件：

```text
chromosome_or_scaffold    gene_id    start    end
```

该列序对应本仓库准备脚本和本文所述的`Quick Run MCScanX Wrapper`路线，不应将ID与染色体列互换。不同MCScanX版本和TBtools-II入口对位置文件格式的约定可能不同；若改用`One Step MCScanX`或直接运行其他版本，应依实际说明核实，不要将标准九列GFF3、BED文件和这里的简化四列文件视为可互换格式。

脚本会打印`Protein records in FASTA`、`Protein IDs with coordinates in GFF`和`Rows written`三项计数。前者是输入蛋白质记录数，第二项是从注释匹配到位置的蛋白质ID数，第三项是实际写入位置文件的行数。`Rows written`必须大于0；如果它远低于输入蛋白记录数，先调查注释中是否缺少蛋白质、ID字段选择错误或异构体策略不同，不要继续运行DIAMOND。无异常重复映射时，后两项通常应一致。生成后马上检查文件是否存在、四列是否有效、ID是否重复：

```bash
ls -lh data/intermediate/TK.clean.protein.fasta data/intermediate/TK.gff
grep -c '^>' data/intermediate/TK.clean.protein.fasta
wc -l data/intermediate/TK.gff
awk -F '\t' 'NF != 4 || $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ {print "格式异常，行号：" NR; bad=1} END {if (!bad) print "四列及坐标格式检查通过"; exit bad}' data/intermediate/TK.gff
cut -f2 data/intermediate/TK.gff | sort | uniq -d
head -n 3 data/intermediate/TK.gff
```

四列格式检查应显示通过；重复ID命令在没有重复时不输出内容。`wc -l`与蛋白质FASTA序列数不必完全相同，因为有些蛋白质可能无法在指定GFF feature中找到坐标，但差异应能被解释。可用以下命令列出未获得位置的蛋白质ID，结果用于排查，不会改写任何输入文件：

```bash
grep '^>' data/intermediate/TK.clean.protein.fasta | sed 's/^>//; s/[[:space:]].*$//' | sort -u > data/intermediate/protein_ids.txt
cut -f2 data/intermediate/TK.gff | sort -u > data/intermediate/position_ids.txt
comm -23 data/intermediate/protein_ids.txt data/intermediate/position_ids.txt | head -n 20
```

如果未匹配ID很多，先处理注释和FASTA的映射问题。还应抽查染色体/scaffold名与坐标是否处于合理范围；本仓库脚本会将同一ID重复出现的CDS区间汇总为最小起点至最大终点，因此必须确认重复记录确实属于同一目标基因或蛋白质，而不是由错误ID字段造成的跨基因合并。

## 五、DIAMOND蛋白质全对全比对

以下首段命令适用于单物种种内比对：DIAMOND的输入必须是上一步生成、且与位置文件ID一致的`TK.clean.protein.fasta`。先为该FASTA建立DIAMOND数据库，再将同一FASTA作为查询执行蛋白质比对：

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

可以逐条粘贴执行，也可以整段执行；如果`diamond makedb`报错，应先停止并修复数据库步骤，不要继续运行后续比对命令。成功时数据库命令会正常结束，并在中间目录生成`TK.dmnd`。

第一条命令生成数据库文件`TK.dmnd`；第二条命令输出DIAMOND兼容的BLAST tabular格式，默认前两列是query和subject ID，后续列为相似度、比对长度、E-value、bit score等信息。运行完应确认输出非空并抽查：

```bash
ls -lh data/intermediate/TK.dmnd data/intermediate/TK.blast
wc -l data/intermediate/TK.blast
head -n 5 data/intermediate/TK.blast
awk -F '\t' 'NF < 12 {print "列数不足，行号：" NR; bad=1} END {if (!bad) print "比对表字段检查通过"; exit bad}' data/intermediate/TK.blast
awk -F '\t' 'FNR == NR {valid[$2]=1; next} !($1 in valid) || !($2 in valid) {print "位置文件中找不到比对ID，行号：" FNR; bad=1} END {if (!bad) print "比对ID与位置文件匹配"; exit bad}' data/intermediate/TK.gff data/intermediate/TK.blast
```

输出文件为空时，检查FASTA是否含有效蛋白质序列、数据库是否成功建立、命令路径是否正确，以及DIAMOND运行日志中的错误。`--evalue 1e-5`和每条查询最多保留5个匹配仅为示例参数，不是普适标准；更严格的E-value或较小的`--max-target-seqs`会减少候选匹配，可能降低共线性识别能力。反过来，保留过多弱匹配也可能引入噪声。应结合物种关系、基因组重复情况、研究目的及MCScanX建议确定参数，并在方法中记录DIAMOND版本、E-value、最大匹配数及输出字段。MCScanX上游说明中的常见建议是为每个查询保留少量高分匹配，但具体阈值仍需依据数据判断，不能机械套用。

### 两物种种间比对的具体做法

先按第4节分别准备两个物种的文件，至少应有`A.clean.protein.fasta`、`A.gff`、`B.clean.protein.fasta`和`B.gff`。确认A、B两套文件各自FASTA与GFF的ID对应正确，再确认合并后不会有重名基因ID或染色体名。随后建立两个独立的DIAMOND数据库：

```bash
diamond makedb --in data/intermediate/A.clean.protein.fasta --db data/intermediate/A
diamond makedb --in data/intermediate/B.clean.protein.fasta --db data/intermediate/B
```

接着分别执行四个方向的搜索：A蛋白对A数据库（A-A种内）、B蛋白对B数据库（B-B种内）、A蛋白对B数据库（A-B种间）、B蛋白对A数据库（B-A种间）。下例保留每条query最多5个target；这是示例上限，不等于严格的双向最佳命中筛选：

```bash
diamond blastp --query data/intermediate/A.clean.protein.fasta --db data/intermediate/A.dmnd --out data/intermediate/A_A.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
diamond blastp --query data/intermediate/B.clean.protein.fasta --db data/intermediate/B.dmnd --out data/intermediate/B_B.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
diamond blastp --query data/intermediate/A.clean.protein.fasta --db data/intermediate/B.dmnd --out data/intermediate/A_B.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
diamond blastp --query data/intermediate/B.clean.protein.fasta --db data/intermediate/A.dmnd --out data/intermediate/B_A.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
```

检查四个文件都已生成并非空，随后按文本顺序合并为单个DIAMOND表格，供MCScanX使用：

```bash
ls -lh data/intermediate/A_A.blast data/intermediate/B_B.blast \
  data/intermediate/A_B.blast data/intermediate/B_A.blast
cat data/intermediate/A_A.blast data/intermediate/B_B.blast \
  data/intermediate/A_B.blast data/intermediate/B_A.blast \
  > data/intermediate/AB.blast
wc -l data/intermediate/AB.blast
head -n 5 data/intermediate/AB.blast
awk -F '\t' 'FNR == NR {valid[$2]=1; next} !($1 in valid) || !($2 in valid) {print "位置文件中找不到比对ID，行号：" FNR; bad=1} END {if (!bad) print "合并比对ID与合并位置文件匹配"; exit bad}' data/intermediate/AB.gff data/intermediate/AB.blast
```

后续在TBtools-II`Quick Run MCScanX Wrapper`中选择合并位置文件`data/intermediate/AB.gff`和合并比对文件`data/intermediate/AB.blast`。如果本次仅做两物种之间的种间区块，而不需要两物种各自的种内区块，可根据MCScanX参数或经验证的筛选流程排除种内命中；不能简单删除文件中的某几列。MCScanX论文建议多物种分析考虑合并各物种自比对与两两比对，并使用研究方案规定的best-hit筛选。DIAMOND的`--max-target-seqs 5`只是每条query保留的target数量上限，不是reciprocal best hit（RBH）；需要RBH时必须额外运行明确的互为最佳匹配筛选并保存筛选规则。

## 六、运行MCScanX并检查共线性结果

本流程使用TBtools-II的`Quick Run MCScanX Wrapper`作为图形界面入口。不同版本的菜单布局或字段名称可能略有变化，因此以工具搜索结果中的完整名称为准，不要把`One Step MCScanX`当作完全相同的入口。先确认`TK.blast`和`TK.gff`均已生成，`TK.blast`是DIAMOND的制表符比对结果而不是`.dmnd`数据库文件，`TK.gff`是本仓库生成的四列位置文件而不是原始九列GFF3。位置文件必须无表头，且两文件中的蛋白质ID完全一致。

在TBtools-II主窗口中按以下顺序操作：

1. 在工具搜索框中输入`Quick Run MCScanX Wrapper`，从搜索结果中打开名称完全一致的工具。若只看到`One Step MCScanX`，不要将其当作同一界面继续填写；先确认TBtools-II版本及工具是否已安装。
2. 在`Input .blast File`栏选择`data/intermediate/TK.blast`。选文件时确认扩展名为`.blast`，不要误选`TK.dmnd`、原始蛋白质FASTA或空文件。
3. 在`Input .gff File`栏选择`data/intermediate/TK.gff`。这里应选择本仓库生成的四列位置文件，而不是`data/raw/genome.gff`或`.gff3`原始注释。
4. 在`Output Directory`栏选择已存在的`results/mcscanx/`文件夹。若文件选择窗口不能输入相对路径，可先在文件管理器中定位到仓库目录，再依次进入`results`和`mcscanx`。确认输入框显示的是目标路径后再开始运行。
5. 点击界面中的运行/开始按钮，等待日志或进度提示完成。运行期间不要关闭TBtools-II；若出现错误，先保存或复制错误信息，不要连续重复启动任务。

正常完成后，输出目录应含以输入前缀`TK`命名的`.collinearity`结果，预期路径为`results/mcscanx/TK.collinearity`。如果该文件没有出现，或输出目录仍为空，先检查界面日志中实际使用的输入路径、文件前缀、权限和输入格式。输入文件由仓库代码生成，并不保证所有TBtools-II版本界面完全相同；字段名称不一致时先对照当前版本帮助说明，不要仅凭位置猜测。

结果目录中应出现以`TK`为前缀的MCScanX结果，核心文件为`TK.collinearity`。先检查输出是否存在、大小是否非零，再查看区块标题和若干基因对：

```bash
ls -lh results/mcscanx
head -n 40 results/mcscanx/TK.collinearity
```

`.collinearity`中的区块标题通常以`## Alignment`开头，后续行包含区块内成对的基因ID。确认这些ID能在`TK.gff`第二列和`TK.blast`前两列中找到，并检查区块是否包含多个连续基因对。输出空、ID大面积不匹配或区块数量明显低于预期时，暂不解释生物学结果；依次回查蛋白质ID清理、注释属性字段、位置文件列序和坐标、DIAMOND筛选参数、目标物种或基因组版本。若直接运行命令行版MCScanX，按该版本的文件命名和参数约定运行，并记录准确命令及版本；不要混用不同软件入口的输入格式。

若做两物种比较，输入栏分别改选`data/intermediate/AB.blast`和`data/intermediate/AB.gff`，输出文件前缀应与合并文件前缀一致，预期核心结果为`results/mcscanx/AB.collinearity`；不要仍选择单物种的`TK`文件。

## 七、筛选目标基因相关的共线性基因对

目标列表中的ID必须与`TK.gff`第二列和`.collinearity`中的ID完全一致。将示例复制到本次运行的中间目录，再编辑副本；每行填写一个目标基因或蛋白质ID，不加列名、逗号或说明文字。空行和以`#`开头的注释行会被忽略：

```bash
cp -n configs/targets.example.txt data/intermediate/targets.txt
```

`-n`表示目标文件已存在时不覆盖；若之前已经创建过`targets.txt`，直接编辑该文件即可。

下面的`GENE000001`等仅为示例占位ID，须替换成你实际要分析的ID；如果原样保留，筛选脚本只会尝试查找这些示例字符串，通常不会命中真实数据。

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
  --targets data/intermediate/targets.txt \
  --out_prefix results/mcscanx/target_collinearity
```

脚本会生成`results/mcscanx/target_collinearity_pairs.tsv`和`results/mcscanx/target_collinearity_summary.tsv`。建议检查文件是否存在并查看表头及前几行：

```bash
ls -lh results/mcscanx/target_collinearity_pairs.tsv results/mcscanx/target_collinearity_summary.tsv
head -n 5 results/mcscanx/target_collinearity_pairs.tsv
head -n 5 results/mcscanx/target_collinearity_summary.tsv
```

`pairs.tsv`逐条记录目标基因、伙伴基因、区块、坐标、关系及区块标题；`summary.tsv`汇总目标基因是否命中、伙伴基因数、区块数及同染色体或异染色体关系。`num_partner_genes`统计唯一伙伴基因数，而`num_pairs`按结果行计数，同一对基因若出现在多个区块可能重复出现。因此，报告“不同基因对数”时不能直接把`num_pairs`当作去重后的基因对数量。

若要生成不区分方向的唯一基因对表，保留原始`pairs.tsv`，并额外运行下列命令。命令按前两列ID的字典顺序建立唯一键，同一对基因即使在不同区块重复出现也只在新文件中保留一次；区块和坐标保留首次出现记录，原始文件不被修改：

```bash
awk -F '\t' 'BEGIN {OFS="\t"} NR==1 {print; next} {
  a=$1; b=$2
  if (a>b) {tmp=a; a=b; b=tmp}
  key=a SUBSEP b
  if (!seen[key]++) print
}' results/mcscanx/target_collinearity_pairs.tsv \
  > results/mcscanx/target_collinearity_unique_pairs.tsv
```

`found_in_collinearity = NO`表示筛选脚本没有在所解析的共线性区块中找到目标ID，不等于目标基因不存在或没有任何同源关系。先确认目标ID是否存在于位置文件，再核对`.collinearity`格式、注释版本和ID层级。该脚本按MCScanX标准`.collinearity`文本结构解析；若使用其他工具或不同输出格式，应先确认原始文件与解析逻辑兼容。

## 八、准备Ka/Ks基因对列表

Ka/Ks针对一对编码序列估算，不是为单个基因计算一个独立常数。每一行须代表研究设计中要比较的一对同源基因或复制基因，并且两端必须对应可比的编码序列ID。候选基因对可来自上一步共线性结果，也可来自已完成的同源分析；应预先说明纳入标准，避免把所有序列任意两两组合。若候选对来自唯一共线性基因对表，可直接提取前两列：

```bash
tail -n +2 results/mcscanx/target_collinearity_unique_pairs.tsv | cut -f1,2 \
  > data/intermediate/kaks_pairs.tsv
```

若候选对不是来自该唯一对表，则将示例复制到`data/intermediate/kaks_pairs.tsv`并编辑。文件必须为纯文本，每行恰好两个ID，由真实Tab字符分隔且不加表头。下文中的`ID_A<Tab>ID_B`仅表示两个ID之间需要一个Tab字符；实际文件中不要输入尖括号或`Tab`三个字母：

```text
ID_A<Tab>ID_B
ID_C<Tab>ID_D
```

检查每行是否恰有两列，再生成候选CDS的唯一ID清单：

```bash
awk -F '\t' 'NF != 2 || $1 == "" || $2 == "" || $1 == $2 {print "列数、空ID或自配对错误，行号：" NR; bad=1} END {if (NR == 0) {print "基因对文件为空"; exit 1} if (!bad) print "基因对文件格式检查通过"; exit bad}' data/intermediate/kaks_pairs.tsv
awk -F '\t' '{print $1; print $2}' data/intermediate/kaks_pairs.tsv | sort -u \
  > data/intermediate/kaks_ids.txt
wc -l data/intermediate/kaks_pairs.tsv data/intermediate/kaks_ids.txt
head -n 5 data/intermediate/kaks_pairs.tsv
```

若格式检查报错，应检查是否误用空格或逗号代替Tab、是否包含表头、空行或多余列。基因对数可能大于唯一ID数的一半，因为同一基因可参与多对比较。确认每对序列的物种、关系、来源区块和筛选依据，并将这些信息与结果关联保存。共线性是基因关系判断的证据之一，并非直系同源判定的充分条件；种内复制基因对与跨物种基因对的演化解释不能混为一类。

## 九、提取候选CDS

若已有全量CDS FASTA，先用前述ID清单提取候选序列。该命令假设FASTA标题的第一个字段与`kaks_ids.txt`完全一致：

```bash
mkdir -p results/kaks
python3 scripts/extract_cds_by_ids.py \
  --cds data/raw/all.cds.fa \
  --ids data/intermediate/kaks_ids.txt \
  --out results/kaks/target.cds.fa
```

若没有CDS FASTA，但有与注释版本一致的基因组FASTA和GFF/GFF3，可先用gffread提取全量CDS，再从中筛选候选序列。`-g`指定基因组FASTA，`-x`指定CDS输出路径，最后一个参数是注释文件：

```bash
gffread data/raw/genome.gff \
  -g data/raw/genome.fasta \
  -x data/intermediate/all.cds.fa

python3 scripts/extract_cds_by_ids.py \
  --cds data/intermediate/all.cds.fa \
  --ids data/intermediate/kaks_ids.txt \
  --out results/kaks/target.cds.fa
```

gffread要求注释中的序列名称能够在基因组FASTA标题中找到，且注释需具有足以还原转录本外显子结构的记录和父子关系。若seqid不一致、注释缺少转录本/CDS关系或参考版本不匹配，提取结果可能为空、不完整或错误。运行后先查看输出大小和FASTA标题。gffread输出ID可能是转录本ID，而候选对使用的可能是基因ID或蛋白质ID；必须先核对映射，不能直接假定相同：

```bash
grep -m 5 '^>' data/intermediate/all.cds.fa
head data/intermediate/kaks_ids.txt
ls -lh data/intermediate/all.cds.fa
```

如果FASTA标题与配对清单看起来相似但不完全相同，可先抽取全量CDS的标题ID，再列出候选清单中找不到的ID：

```bash
grep '^>' data/intermediate/all.cds.fa | sed 's/^>//; s/[[:space:]].*$//' | sort -u \
  > data/intermediate/all_cds_ids.txt
comm -23 data/intermediate/kaks_ids.txt data/intermediate/all_cds_ids.txt | head -n 20
```

`comm`没有输出表示所有请求ID都出现在CDS FASTA中；若有输出，说明这些ID不能被精确匹配。先查阅注释的gene-transcript-protein对应关系，明确转换ID的方法，再决定是否生成统一的配对/序列ID，不要靠截短ID来消除差异。

提取脚本会报告找到和缺失的ID，并写出`results/kaks/target.cds.fa`。`Found`应与请求ID数量相符，且`Missing`应为空；若存在缺失，先回查候选ID列表、CDS FASTA标题、版本号和代表转录本映射。不要仅通过模糊字符串替换强行对应。还应检查提取条目数：

```bash
grep -c '^>' results/kaks/target.cds.fa
wc -l data/intermediate/kaks_ids.txt
grep '^>' results/kaks/target.cds.fa | sed 's/^>//; s/[[:space:]].*$//' | sort -u \
  > data/intermediate/extracted_cds_ids.txt
comm -23 data/intermediate/kaks_ids.txt data/intermediate/extracted_cds_ids.txt | head -n 20
```

前两个计数通常应相等；若CDS FASTA内同一ID重复出现，提取条目数可能异常，应先确认输入FASTA记录唯一性。最后的`comm`命令列出“请求但未提取到”的ID；正确完成时应没有输出。若仍有ID出现，暂不运行Ka/Ks，回查提取脚本的`Found/Missing`提示和ID映射。

## 十、CDS质量检查

对候选CDS运行仓库提供的基础检查，并将检查表保存到结果目录：

```bash
python3 scripts/validate_cds_lengths.py \
  --cds results/kaks/target.cds.fa \
  > results/kaks/cds_validation.tsv
head -n 10 results/kaks/cds_validation.tsv
```

检查表各列为序列ID、长度、长度是否为3的倍数、是否以ATG起始、是否以终止密码子结束，以及是否含ATCGN以外字符。该脚本只进行基础格式和长度提示，不会自动翻译序列，也不会检测所有内部终止密码子、移码、错误剪接、外显子边界或密码子比对问题。含N的序列不一定不合法，但模糊碱基会降低替换率估计的确定性；不同遗传密码表、叶绿体/线粒体基因及部分不完整CDS也可能不符合标准ATG起始或终止规则。因此，该表不是自动筛除规则，而是逐条复核的检查清单。

对每一对待分析序列至少确认ID正确对应同源编码序列；CDS长度合理且阅读框完整或其不完整状态已知；没有未解释的内部终止或移码；双方使用一致的遗传密码表；序列比对按密码子方式建立，缺口和低质量区域已妥善处理。若序列高度相似或差异很大，也要检查是否选错异构体、物种或同源伙伴。完成必要复核后再进入Ka/Ks计算。

## 十一、计算Ka/Ks

本仓库不包含Ka/Ks估算程序。以下以TBtools-II中的`Simple Ka/Ks Calculator (NG)`为例；不同版本的界面和输出字段可能变化，应以实际安装版本的说明为准。计算前准备两份输入：`results/kaks/target.cds.fa`中必须有候选基因对两端所需的全部CDS，`data/intermediate/kaks_pairs.tsv`中每行必须恰有两个ID且无表头。两个文件里的ID必须完全一致：配对表第一列和第二列分别指向FASTA标题`>`后第一个空格之前的ID。`Input GenePair File`必须选择这个两列配对清单，不能选择原始`.collinearity`文件，也不能直接选择带表头和多列的`target_collinearity_pairs.tsv`。

若手动编辑配对表，用纯文本编辑器新建或打开`data/intermediate/kaks_pairs.tsv`，输入第一条ID，按一次键盘Tab，再输入第二条ID；按Enter换到下一行，逐对填写。不要输入列名，不要用空格或逗号替代Tab。若候选对来自上一节生成的唯一共线性基因对表，优先使用前面给出的`tail`和`cut`命令生成，避免手工录入错误。运行以下命令可检查Tab是否真实存在：

```bash
sed -n '1,5l' data/intermediate/kaks_pairs.tsv
```

正确的两列通常会显示为`ID_A\tID_B$`，其中`\t`是Tab的可见表示，`$`表示该行结束。若显示`ID_A ID_B$`且中间是空格，或显示文字`<Tab>`，应返回文本编辑器修正。

在TBtools-II中按以下顺序操作：

1. 在工具搜索框中输入`Simple Ka/Ks Calculator (NG)`，打开名称完全一致的计算器。确认本次确实选择了NG计算方法；若研究方案要求其他方法，应改选对应工具并在方法中另行记录。
2. 在`Input CDS File`栏选择`results/kaks/target.cds.fa`。这是已筛选的候选CDS文件，不是基因组FASTA，也不是蛋白质FASTA。
3. 在`Input GenePair File`栏选择`data/intermediate/kaks_pairs.tsv`。确认它是无表头的两列Tab分隔文件。
4. 在`Output Table File`栏指定`results/kaks/kaks.result.tsv`。若保存对话框只要求选择文件夹，则先选择`results/kaks/`，再按界面要求填写输出文件名。
5. 对照三个输入框显示的路径逐项检查后，点击运行/开始。等待任务完成；若工具报告ID缺失、序列缺失或格式不符，先修复输入文件再重跑，不要把未输出的基因对当作Ka/Ks为0。

计算完成后应看到`results/kaks/kaks.result.tsv`。运行下面的命令确认文件非空，并检查列名、前几条记录和是否出现NaN/Inf：

```bash
ls -lh results/kaks/kaks.result.tsv
head -n 10 results/kaks/kaks.result.tsv
wc -l results/kaks/kaks.result.tsv data/intermediate/kaks_pairs.tsv
grep -Ein 'nan|inf' results/kaks/kaks.result.tsv || true
```

逐行确认结果中的配对ID与输入一致，并检查空值或明显异常值。结果行数通常应与输入配对数相近；若少于输入配对数，检查哪些配对被跳过以及软件日志，不能把缺行解释为Ka/Ks为0。保存TBtools-II版本、计算器名称、模型、遗传密码表、参数、序列对齐文件（若工具可输出）及运行日期。NG通常指Nei–Gojobori类估算方法；不同比对策略和替代模型可能产生不同结果。未经检查的CDS直接逐碱基比较可能因插入缺失、移码、错误转录本或注释不完整而产生无效估计，因此应确认比对适合密码子层面分析。

## 十二、结果判读与报告

| Ka/Ks或ω | 常见初步解释 | 注意事项 |
|---|---|---|
| `< 1` | 与纯化选择相符 | 需排除序列错配、错误CDS及Ks估算问题 |
| `≈ 1` | 与中性演化相符 | 仅是简化解释，受模型和统计不确定性影响 |
| `> 1` | 可能提示正选择或功能分化 | 单个基因对不足以证明适应性进化，应结合序列、模型和多基因证据 |
| `NaN`或`Inf` | 当前数据或模型下无法稳定估算 | 不代表纯化选择，也不应直接替换为0 |

建议按下表整理最终结果，分别报告基因数和基因对数，不要把二者混为一谈：

| 指标 | 推荐口径 | 检查或计算方式 |
|---|---|---|
| 输入蛋白质数 | 蛋白质FASTA中的唯一ID数 | 以FASTA标题首字段去重统计；必要时另报原始条目数 |
| 有效位置ID数 | `TK.gff`第二列中的唯一ID数 | 从第二列提取ID并去重计数 |
| 共线性基因数 | 至少出现在一个共线性区块中的唯一基因ID数 | 从`.collinearity`基因对行提取两端ID后合并去重 |
| 共线性基因比例 | 共线性唯一基因数÷进入MCScanX的位置文件唯一ID数 | 明确使用位置文件ID作分母；若分母用蛋白质总数，需另行注明覆盖率 |
| 目标基因命中数 | `found_in_collinearity`为YES的不同目标基因数 | 以目标汇总表为准，并核对实际字段和值 |
| 目标相关唯一基因对数 | 按无向基因对去重后的行数，不含表头 | 使用`target_collinearity_unique_pairs.tsv`，不要使用原始重复行数 |
| Ka/Ks有效结果数 | 同时具有可解释Ka、Ks和Ka/Ks的基因对数 | 统计时将NaN、Inf及其他失败结果单独列出，不要当作0 |

如果当前MCScanX输出的基因对行遵循标准`.collinearity`格式，可用以下命令从区块配对行提取并统计唯一共线性基因ID；运行后抽查文件内容，确认提取到的是基因ID而非索引号：

```bash
grep '^>' data/intermediate/TK.clean.protein.fasta | sed 's/^>//; s/[[:space:]].*$//' | sort -u | wc -l
cut -f2 data/intermediate/TK.gff | sort -u | wc -l
awk '/^[[:space:]]*[0-9]+-[[:space:]]+[0-9]+:/ {print $3; print $4}' \
  results/mcscanx/TK.collinearity | sort -u > data/intermediate/collinear_gene_ids.txt
wc -l data/intermediate/collinear_gene_ids.txt
head data/intermediate/collinear_gene_ids.txt
```

以上命令依次输出唯一蛋白质ID数、唯一位置ID数，并生成唯一共线性基因ID清单。将共线性基因数除以位置ID数，即为上述口径下的共线性基因比例。若MCScanX输出格式不同，先根据实际区块配对行调整提取方式，避免直接套用导致把索引当ID。对跨物种分析，应分别报告每个物种中进入位置文件和出现在共线性区块的基因数，不能用合并总数掩盖某一物种的数据缺失。

Ka/Ks是针对一对同源编码序列的估计值，不是单个基因固有属性。Ks较高可能发生替换饱和；Ka或Ks接近0时比值可能不稳定。报告每一对的原始Ka、Ks和Ka/Ks，并标记无法估算的记录。若`Ka=0`且`Ks=0`，Ka/Ks为`0/0`，在数学上未定义，不应替换成0。Ka/Ks小于1、大于1或接近1都只是初步解释，不能仅凭单个比值下结论。

方法记录至少应包括：比较类型（种内或种间）、物种及参考基因组/注释版本、蛋白质与CDS来源、ID映射和代表转录本规则、DIAMOND/MCScanX/TBtools-II/gffread版本、DIAMOND参数、输入及匹配数量、共线性筛选规则、CDS质量检查方法、Ka/Ks计算器/模型/密码表/参数、异常结果处理方式，以及结果统计口径。保留原始输入、关键中间文件和运行记录，可使后续复核不必从头猜测使用过的参数。

各阶段预期文件可按此核对：

| 阶段 | 主要文件 | 继续下一步前的最低检查 |
|---|---|---|
| 输入准备 | `data/raw/`中的蛋白质FASTA、GFF/GFF3及CDS或基因组FASTA | 文件非空、版本一致、ID层级已确认 |
| MCScanX输入生成 | `TK.clean.protein.fasta`、`TK.gff` | ID匹配情况可解释，位置文件四列且ID不重复 |
| DIAMOND | `TK.dmnd`、`TK.blast` | 数据库和比对文件非空，比对表字段及ID正常 |
| MCScanX | `TK.collinearity` | 文件非空，抽查区块基因ID与输入一致 |
| 目标筛选 | `target_collinearity_pairs.tsv`、`target_collinearity_summary.tsv` | 目标命中状态、伙伴基因和区块数量合理 |
| CDS准备 | `target.cds.fa`、`cds_validation.tsv` | 候选对所需ID均找到，异常CDS已复核 |
| Ka/Ks | `kaks.result.tsv` | 输入配对与输出配对一致，异常值有记录 |

## 十三、常见问题

**`conda`、`diamond`或`gffread`提示找不到命令。** 先确认终端已激活环境`collinearity`，再用`conda list`检查软件是否安装；若刚安装完仍找不到，关闭并重新打开终端后重新激活环境。不要在未激活环境时误判为软件未安装。

**`Rows written`为0。** 通常表示`--feature_type`或`--protein_attr`与GFF实际内容不一致。检查GFF第三列feature类型和第九列属性字段，并确认提取的ID与蛋白质FASTA标题一致。尤其注意脚本默认使用`CDS`和`Protein_Accession`，且属性解析按`key=value`格式进行。

**只有少量蛋白质写入`TK.gff`。** 检查蛋白质FASTA与注释中的ID是否只部分匹配，重点核对版本号、前缀、转录本后缀、空格后的描述及ID层级。

**MCScanX只识别出很少基因。** 检查位置文件是否为Tab分隔的四列，且按本流程约定的“染色体或scaffold、gene ID、start、end”顺序排列；不要将gene ID放在第一列、染色体放在第二列。另须确认第三、四列为有效数值坐标、每个基因ID只对应一个位置、ID与比对文件匹配，且输入蛋白质及匹配阈值符合研究设计。

**DIAMOND比对文件为空或过小。** 检查FASTA记录数及序列内容、数据库文件是否生成、查询文件和数据库是否来自同一批ID；再确认E-value和最大匹配数设置是否过严。跨物种分析时检查两物种蛋白质FASTA是否都纳入合并输入，且ID唯一。参数调整应记录理由，不要只为增加区块数而无依据地放宽阈值。

**MCScanX没有生成`.collinearity`或结果目录为空。** 查看TBtools-II运行日志和输出目录，确认比对文件、四列位置文件和文件前缀与界面要求相符。若通过GUI运行后仍无结果，先用少量测试数据验证入口和格式，再用完整数据；不要反复点击运行并覆盖之前的日志。

**TBtools-II或MCScanX提示ID不一致。** 确认`TK.clean.protein.fasta`标题、`TK.blast`前两列和`TK.gff`第二列采用兼容且唯一的ID。

**TBtools-II的`One Step MCScanX`报告ID不一致。** 可将流程拆开定位问题：先检查DIAMOND生成的比对表格，再用本仓库脚本准备四列位置文件，最后将两份文件分别交给`Quick Run MCScanX Wrapper`。这是一种排错路径，并不表示所有版本的`One Step MCScanX`都不兼容；切换入口时应核实其文件格式和ID字段要求。

**目标基因均未命中共线性结果。** 确认目标列表使用的ID层级与`TK.gff`及`.collinearity`完全一致，并确认目标ID存在于输入数据中。

**CDS提取出现Missing。** 回到注释和FASTA标题建立明确的gene、transcript、protein映射，并检查是否使用了正确的代表转录本。

**gffread输出为空或CDS数异常少。** 核对GFF第一列seqid是否与基因组FASTA标题一致、GFF与FASTA是否来自同一组装版本、注释是否含有正确的转录本/CDS父子关系。gffread提取的CDS ID与配对文件可能不是同一层级，应先处理ID映射，不要把空结果交给Ka/Ks计算器。

**Ka/Ks结果为NaN或Inf。** 若`Ka=0`且`Ks=0`，比值为`0/0`，在数学上未定义；该基因对不能据此进行选择压力判断，也不应将结果改写为0。其他NaN或Inf还可能与有效替换位点不足、序列相同或近似相同、单独的Ka或Ks为0、阅读框错误、密码子比对异常或模型不适用有关。应逐对检查原始CDS、密码子比对和模型，并在结果表中记录异常原因。

## 十四、数据安全与提交前检查

`.gitignore`默认排除FASTA、GFF/GTF、GenBank文件、DIAMOND数据库与比对结果、MCScanX结果、常见表格、图片、日志和报告，以及`data/raw/`、`data/intermediate/`和`results/`中的实际内容。忽略规则不能替代提交前审查；上传前运行：

```bash
git status
```

确认没有真实序列、未公开数据、样本信息或不宜公开的分析结果被加入提交。

正式汇报前还应完成以下复核：输入数据版本及来源已记录；所有工具版本和参数可追溯；ID映射及代表转录本规则一致；共线性输入与输出的计数口径明确；目标基因对经过去重并保留原始筛选表；CDS异常、NaN、Inf及Ks可能饱和的结果均已单独标记；统计结论没有把共线性等同于直系同源，也没有把单个Ka/Ks大于1直接表述为正选择的证明。

## 十五、Windows原生环境（PowerShell）操作路线

本节面向Windows11原生PowerShell，不使用Bash语法。分析原理、输入要求、ID核对、结果解释及参数选择仍遵循前文；本节补充Windows中的安装方式、命令写法和文件检查方法。若软件或工具封装在当前Windows版本不可用，应暂停对应步骤，改用第2至第14节所述Ubuntu/Linux环境或WSL，而不是把Linux命令直接粘贴到PowerShell。

### 15.1 Windows路线的适用范围与软件准备

Windows原生路线可用PowerShell运行本仓库Python脚本，并通过DIAMOND官方Windows独立程序完成蛋白质比对；TBtools-II可安装Windows版本并用于图形界面操作。需要注意两项边界：MCScanX上游命令行说明以Linux和macOS为运行环境，因此Windows路线应以TBtools-II Windows版的`Quick Run MCScanX Wrapper`确实存在且成功运行为前提；如果该入口缺失、报错或生成结果异常，应改用Linux/WSL运行MCScanX。gffread的Bioconda平台列表不包含Windows；若只有基因组FASTA和GFF/GFF3、尚无CDS FASTA，最稳妥的处理是使用Linux/WSL路线提取CDS，再按同一ID规则继续分析。不要在Windows上照抄`conda create ... gffread mcscanx`并假设这些包一定可用。

先打开Windows PowerShell，检查Python启动器和Git：

```powershell
py -0p
py -3 --version
git --version
```

若`py -3`不可用但`python --version`能显示Python3，则本节的`py -3`命令可改为`python`。本仓库脚本只依赖Python标准库，不需要安装额外Python包。建议从DIAMOND[官方发布页](https://github.com/bbuchfink/diamond/releases)下载Windows压缩包，将其中的`diamond.exe`放到项目的`tools\diamond\diamond.exe`。DIAMOND不是通过`pip install`安装；若启动时报缺少`VCRUNTIME`或类似运行库错误，应从Microsoft官方渠道安装适用于系统的Visual C++ Redistributable，然后重新打开PowerShell验证。不要双击`diamond.exe`启动，它需要在终端中附带参数运行。

从[TBtools-II官方发布页](https://github.com/CJ-Chen/TBtools-II/releases)下载Windows安装包并完成安装。打开软件后，确认主界面可启动，并在工具搜索框中检查是否存在`Quick Run MCScanX Wrapper`及`Simple Ka/Ks Calculator (NG)`。记录软件版本；若Windows版未提供可用的MCScanX封装，不要把“TBtools已安装”当作MCScanX已可运行。

首次下载仓库时运行`git clone`，随后进入项目根目录；若仓库已经下载，只执行`Set-Location`并将示例路径替换为本机实际位置：

```powershell
git clone https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline.git
Set-Location 'D:\Bioinfo-collinearity-kaks-pipeline'
```

上例中的`D:\Bioinfo-collinearity-kaks-pipeline`是路径示例，必须替换成实际克隆目录。随后创建工作目录并检查当前位置：

```powershell
New-Item -ItemType Directory -Force -Path .\data\raw, .\data\intermediate, .\results\mcscanx, .\results\kaks, .\tools\diamond | Out-Null
Get-Location
```

将输入文件放入`data\raw\`，将官方压缩包中的`diamond.exe`放入`tools\diamond\`。该可执行文件仅供本机调用，提交Git前确认没有将其加入暂存区。路径中尽量避免中文、空格及特殊符号；若项目路径含空格，在命令参数中用引号包住路径。PowerShell命令均假定从项目根目录执行。

### 15.2 检查原始文件及ID

用以下PowerShell命令确认主要输入文件存在，并查看蛋白质FASTA标题及GFF/GFF3记录：

```powershell
Get-Item .\data\raw\protein.fasta, .\data\raw\genome.gff | Select-Object FullName, Length
Select-String -Path .\data\raw\protein.fasta -Pattern '^>' | Select-Object -First 5
Select-String -Path .\data\raw\genome.gff -Pattern 'Protein_Accession=' | Select-Object -First 5
```

最后一条只是按文本查找属性示例；仍须人工确认命中的GFF记录第三列确为预期feature（例如`CDS`），并确认属性值与FASTA标题首字段完全一致。如果真实属性是`protein_id=`，将搜索字符串和后续脚本参数改为该字段。查看基因组标题时可运行`Select-String -Path .\data\raw\genome.fasta -Pattern '^>' | Select-Object -First 5`；该文件仅在后续需要gffread时使用。

### 15.3 生成MCScanX输入

在PowerShell中，Python脚本路径使用Windows相对路径，并以`py -3`启动。默认属性字段为`Protein_Accession`时运行：

```powershell
py -3 .\scripts\prepare_mcscanx_inputs.py --gff .\data\raw\genome.gff --protein .\data\raw\protein.fasta --prefix TK --outdir .\data\intermediate
```

若蛋白质ID位于GFF的`protein_id`属性，改为：

```powershell
py -3 .\scripts\prepare_mcscanx_inputs.py --gff .\data\raw\genome.gff --protein .\data\raw\protein.fasta --prefix TK --outdir .\data\intermediate --feature_type CDS --protein_attr protein_id
```

若实际ID位于`mRNA`记录的`ID`属性，可在确认其确与蛋白质FASTA ID一致且坐标代表目标转录本后，使用`--feature_type mRNA --protein_attr ID`。`--prefix TK`只影响输出文件名，不会修改序列ID。种间分析仍需分别处理A、B两套输入，确保蛋白质ID及染色体/scaffold名在合并前已全局唯一：

```powershell
py -3 .\scripts\prepare_mcscanx_inputs.py --gff .\data\raw\A.gff --protein .\data\raw\A.protein.fasta --prefix A --outdir .\data\intermediate
py -3 .\scripts\prepare_mcscanx_inputs.py --gff .\data\raw\B.gff --protein .\data\raw\B.protein.fasta --prefix B --outdir .\data\intermediate
```

单物种运行后，检查输出存在、大小非零、FASTA记录数及位置文件前几行：

```powershell
Get-Item .\data\intermediate\TK.clean.protein.fasta, .\data\intermediate\TK.gff | Select-Object FullName, Length
(Select-String -Path .\data\intermediate\TK.clean.protein.fasta -Pattern '^>').Count
(Get-Content .\data\intermediate\TK.gff | Measure-Object -Line).Lines
Get-Content .\data\intermediate\TK.gff -TotalCount 3
```

位置文件应为无表头的四列Tab分隔文本，列序是染色体/scaffold、ID、起始坐标、终止坐标。可逐行检查字段数和坐标是否为正整数：

```powershell
$lineNumber = 0
$badLine = $false
Get-Content .\data\intermediate\TK.gff | ForEach-Object {
    $lineNumber++
    $fields = $_.Split([char[]]@([char]9))
    if ($fields.Length -ne 4 -or $fields[2] -notmatch '^\d+$' -or $fields[3] -notmatch '^\d+$') {
        Write-Error "位置文件格式异常，行号：$lineNumber"
        $badLine = $true
    }
}
if ($badLine) { throw '请修复位置文件后再继续。' } else { '四列及坐标格式检查通过。' }
```

脚本日志中的`Rows written`须大于0；若明显少于输入蛋白质记录数，应先核对GFF feature类型、属性键、ID层级及转录本规则。重复ID和未匹配ID也应先调查，不能仅凭输出文件存在就进入比对步骤。若输出行数异常，参照第4节的ID核对说明回查原始文件。

还可检查位置文件中是否有重复ID。正确情况下不应出现重复位置ID；若有重复，先回查注释属性键和异构体处理规则：

```powershell
$seenPositionIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
$duplicateIds = [System.Collections.Generic.List[string]]::new()
foreach ($line in [System.IO.File]::ReadLines((Resolve-Path .\data\intermediate\TK.gff).Path)) {
    $fields = $line.Split([char[]]@([char]9))
    if ($fields.Length -ge 2 -and -not $seenPositionIds.Add($fields[1])) { $duplicateIds.Add($fields[1]) }
}
if ($duplicateIds.Count -eq 0) { '未发现重复位置ID。' } else { $duplicateIds | Select-Object -First 20 }
```

### 15.4 使用DIAMOND进行蛋白质比对

先确认DIAMOND可执行文件在预期位置，并查看版本：

```powershell
& .\tools\diamond\diamond.exe version
```

PowerShell中，调用当前目录或相对路径中的`.exe`时使用调用运算符`&`。单物种种内比对依次建立数据库并运行`blastp`：

```powershell
& .\tools\diamond\diamond.exe makedb --in .\data\intermediate\TK.clean.protein.fasta --db .\data\intermediate\TK
if ($LASTEXITCODE -ne 0) { throw 'DIAMOND makedb运行失败。' }
& .\tools\diamond\diamond.exe blastp --db .\data\intermediate\TK.dmnd --query .\data\intermediate\TK.clean.protein.fasta --out .\data\intermediate\TK.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
if ($LASTEXITCODE -ne 0) { throw 'DIAMOND blastp运行失败。' }
```

检查数据库和比对文件是否存在、大小是否非零，并查看前几条匹配：

```powershell
Get-Item .\data\intermediate\TK.dmnd, .\data\intermediate\TK.blast | Select-Object FullName, Length
Get-Content .\data\intermediate\TK.blast -TotalCount 5
```

DIAMOND默认`outfmt 6`应为12列；需要机器检查时运行：

```powershell
$badLine = $false
$lineNumber = 0
Get-Content .\data\intermediate\TK.blast | ForEach-Object {
    $lineNumber++
    if ($_.Split([char[]]@([char]9)).Length -lt 12) {
        Write-Error "比对表列数不足，行号：$lineNumber"
        $badLine = $true
    }
}
if ($badLine) { throw '请检查DIAMOND输出格式。' } else { '比对表字段检查通过。' }
```

再确认DIAMOND比对两端的ID都存在于位置文件第二列中：

```powershell
$validIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
foreach ($line in [System.IO.File]::ReadLines((Resolve-Path .\data\intermediate\TK.gff).Path)) {
    $fields = $line.Split([char[]]@([char]9))
    if ($fields.Length -ge 2) { [void]$validIds.Add($fields[1]) }
}
$lineNumber = 0
foreach ($line in [System.IO.File]::ReadLines((Resolve-Path .\data\intermediate\TK.blast).Path)) {
    $lineNumber++
    $fields = $line.Split([char[]]@([char]9))
    if ($fields.Length -lt 2 -or -not $validIds.Contains($fields[0]) -or -not $validIds.Contains($fields[1])) { throw "比对文件存在位置文件中未找到的ID，行号：$lineNumber" }
}
'比对ID与位置文件匹配。'
```

`--evalue 1e-5`和`--max-target-seqs 5`仍只是示例参数，须根据物种关系、基因组重复情况和研究问题确定，并记录版本与参数。若做两物种比较，先分别建立A、B数据库，再执行A-A、B-B、A-B和B-A四个方向的搜索：

```powershell
& .\tools\diamond\diamond.exe makedb --in .\data\intermediate\A.clean.protein.fasta --db .\data\intermediate\A
& .\tools\diamond\diamond.exe makedb --in .\data\intermediate\B.clean.protein.fasta --db .\data\intermediate\B
& .\tools\diamond\diamond.exe blastp --query .\data\intermediate\A.clean.protein.fasta --db .\data\intermediate\A.dmnd --out .\data\intermediate\A_A.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
& .\tools\diamond\diamond.exe blastp --query .\data\intermediate\B.clean.protein.fasta --db .\data\intermediate\B.dmnd --out .\data\intermediate\B_B.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
& .\tools\diamond\diamond.exe blastp --query .\data\intermediate\A.clean.protein.fasta --db .\data\intermediate\B.dmnd --out .\data\intermediate\A_B.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
& .\tools\diamond\diamond.exe blastp --query .\data\intermediate\B.clean.protein.fasta --db .\data\intermediate\A.dmnd --out .\data\intermediate\B_A.blast --evalue 1e-5 --max-target-seqs 5 --outfmt 6
if ($LASTEXITCODE -ne 0) { throw '请查看最近一次DIAMOND命令的退出状态及错误信息。' }
```

四个命令须逐一确认成功；上面末尾的退出状态检查只能反映最后执行的一条命令，因此还应检查四个`.blast`文件均存在且非空：

```powershell
Get-Item .\data\intermediate\A_A.blast, .\data\intermediate\B_B.blast, .\data\intermediate\A_B.blast, .\data\intermediate\B_A.blast | Select-Object FullName, Length
```

需要合并时，可用.NET文本流按顺序写入且不添加UTF-8 BOM，避免PowerShell 5.1默认文本编码影响第一条序列ID：

```powershell
$blastFiles = @('.\data\intermediate\A_A.blast', '.\data\intermediate\B_B.blast', '.\data\intermediate\A_B.blast', '.\data\intermediate\B_A.blast')
$outputPath = [System.IO.Path]::GetFullPath('.\data\intermediate\AB.blast')
$writer = [System.IO.StreamWriter]::new($outputPath, $false, [System.Text.UTF8Encoding]::new($false))
try {
    foreach ($file in $blastFiles) {
        $reader = [System.IO.StreamReader]::new([System.IO.Path]::GetFullPath($file))
        try { while ($null -ne ($line = $reader.ReadLine())) { $writer.WriteLine($line) } }
        finally { $reader.Dispose() }
    }
}
finally { $writer.Dispose() }
Get-Item .\data\intermediate\AB.blast | Select-Object FullName, Length
```

合并位置文件`AB.gff`时，同样只能连接列顺序一致且ID全局唯一的A、B文件。可用以下代码按文本流重建合并文件，不会把旧内容附加到新结果：

```powershell
$positionFiles = @('.\data\intermediate\A.gff', '.\data\intermediate\B.gff')
$outputPath = [System.IO.Path]::GetFullPath('.\data\intermediate\AB.gff')
$writer = [System.IO.StreamWriter]::new($outputPath, $false, [System.Text.UTF8Encoding]::new($false))
try {
    foreach ($file in $positionFiles) {
        $reader = [System.IO.StreamReader]::new([System.IO.Path]::GetFullPath($file))
        try { while ($null -ne ($line = $reader.ReadLine())) { $writer.WriteLine($line) } }
        finally { $reader.Dispose() }
    }
}
finally { $writer.Dispose() }
Get-Item .\data\intermediate\AB.gff | Select-Object FullName, Length
```

确认A、B位置文件中的基因ID及染色体/scaffold名称在合并前已全局唯一，再将`AB.blast`和`AB.gff`交给TBtools-II。DIAMOND的`--max-target-seqs 5`不是互为最佳匹配（RBH）筛选。

### 15.5 使用TBtools-II识别共线性区块

本节图形界面步骤与第6节相同，区别是从Windows文件选择器中选取本机路径，例如`D:\Bioinfo-collinearity-kaks-pipeline\data\intermediate\TK.blast`和`D:\Bioinfo-collinearity-kaks-pipeline\data\intermediate\TK.gff`。运行前确认比对表为DIAMOND制表符结果，位置文件为四列文件而非原始九列GFF3，并将输出目录设为项目中的`results\mcscanx\`。工具名称应为`Quick Run MCScanX Wrapper`；若Windows版TBtools-II没有该入口，或运行未产生有效`.collinearity`结果，应改在Linux/WSL执行MCScanX，不要将`.blast`或位置文件改名后假设格式已兼容。种间分析使用时，后续筛选步骤中的`TK`文件名前缀也须替换为`AB`，并指向合并后的文件。

运行后检查结果：

```powershell
Get-Item .\results\mcscanx\TK.collinearity | Select-Object FullName, Length
Get-Content .\results\mcscanx\TK.collinearity -TotalCount 40
```

文件应非空，区块标题通常以`## Alignment`开头，区块配对行中的ID应能在位置文件及比对文件中找到。两物种分析时，TBtools输入改为`AB.blast`和`AB.gff`，输出结果应与合并前缀对应。

### 15.6 筛选目标基因及整理候选基因对

将示例目标列表复制为本次运行文件；若文件已存在则直接编辑，避免覆盖已有目标：

```powershell
if (-not (Test-Path .\data\intermediate\targets.txt)) { Copy-Item .\configs\targets.example.txt .\data\intermediate\targets.txt }
notepad .\data\intermediate\targets.txt
```

在记事本中每行保留一个真实目标ID，不加表头、逗号或说明文字。然后运行筛选脚本：

```powershell
py -3 .\scripts\filter_collinearity_targets.py --collinearity .\results\mcscanx\TK.collinearity --gff .\data\intermediate\TK.gff --targets .\data\intermediate\targets.txt --out_prefix .\results\mcscanx\target_collinearity
Get-Item .\results\mcscanx\target_collinearity_pairs.tsv, .\results\mcscanx\target_collinearity_summary.tsv | Select-Object FullName, Length
Get-Content .\results\mcscanx\target_collinearity_pairs.tsv -TotalCount 5
```

若需按无向基因对去重并保留原始表格列，可运行以下PowerShell代码。该代码将两端ID按区分大小写的序数顺序建立唯一键，保留首次出现的整行记录，原始文件不被改写：

```powershell
$inputPath = [System.IO.Path]::GetFullPath('.\results\mcscanx\target_collinearity_pairs.tsv')
$outputPath = [System.IO.Path]::GetFullPath('.\results\mcscanx\target_collinearity_unique_pairs.tsv')
$reader = [System.IO.StreamReader]::new($inputPath)
$writer = [System.IO.StreamWriter]::new($outputPath, $false, [System.Text.UTF8Encoding]::new($false))
$seenPairs = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
try {
    $header = $reader.ReadLine()
    if ($null -eq $header) { throw '目标基因对文件为空。' }
    $writer.WriteLine($header)
    while ($null -ne ($line = $reader.ReadLine())) {
        $fields = $line.Split([char[]]@([char]9))
        if ($fields.Length -lt 2) { continue }
        $first = $fields[0]; $second = $fields[1]
        if ([string]::CompareOrdinal($first, $second) -gt 0) { $temporary = $first; $first = $second; $second = $temporary }
        $key = $first + [char]31 + $second
        if ($seenPairs.Add($key)) { $writer.WriteLine($line) }
    }
}
finally { $reader.Dispose(); $writer.Dispose() }
Get-Item $outputPath | Select-Object FullName, Length
```

若候选基因对来自该唯一表，可将前两列写成无表头的`kaks_pairs.tsv`。下面的代码假定表头确为`target_gene`和`partner_gene`，脚本会拒绝缺少这两列的输入：

```powershell
$pairTable = '.\results\mcscanx\target_collinearity_unique_pairs.tsv'
$pairRows = Import-Csv -Delimiter "`t" -Path $pairTable
if ($pairRows.Count -eq 0 -or $null -eq $pairRows[0].target_gene -or $null -eq $pairRows[0].partner_gene) { throw '未找到target_gene或partner_gene列，请检查表头。' }
$pairPath = [System.IO.Path]::GetFullPath('.\data\intermediate\kaks_pairs.tsv')
$pairWriter = [System.IO.StreamWriter]::new($pairPath, $false, [System.Text.UTF8Encoding]::new($false))
try { foreach ($row in $pairRows) { $pairWriter.WriteLine("$($row.target_gene)`t$($row.partner_gene)") } }
finally { $pairWriter.Dispose() }
```

如果手动制作`kaks_pairs.tsv`，使用纯文本编辑器；每行两端ID之间必须按一次Tab键，不能有表头。随后逐行核对两列、空ID及自配对，并生成区分大小写的唯一CDS ID列表：

```powershell
$pairPath = [System.IO.Path]::GetFullPath('.\data\intermediate\kaks_pairs.tsv')
$idsPath = [System.IO.Path]::GetFullPath('.\data\intermediate\kaks_ids.txt')
$seenIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
$idsWriter = [System.IO.StreamWriter]::new($idsPath, $false, [System.Text.UTF8Encoding]::new($false))
$lineNumber = 0; $badPair = $false
try {
    foreach ($line in [System.IO.File]::ReadLines($pairPath)) {
        $lineNumber++
        $fields = $line.Split([char[]]@([char]9))
        if ($fields.Length -ne 2 -or [string]::IsNullOrWhiteSpace($fields[0]) -or [string]::IsNullOrWhiteSpace($fields[1]) -or $fields[0] -ceq $fields[1]) {
            Write-Warning "基因对格式异常，行号：$lineNumber"
            $badPair = $true
            continue
        }
        foreach ($id in $fields) { if ($seenIds.Add($id)) { $idsWriter.WriteLine($id) } }
    }
}
finally { $idsWriter.Dispose() }
if ($lineNumber -eq 0) { throw '基因对文件为空。' }
if ($badPair) { throw '请修复基因对文件后再继续。' }
Get-Item $pairPath, $idsPath | Select-Object FullName, Length
```

### 15.7 提取并检查CDS

若已有全量CDS FASTA，运行仓库脚本提取候选CDS：

```powershell
py -3 .\scripts\extract_cds_by_ids.py --cds .\data\raw\all.cds.fa --ids .\data\intermediate\kaks_ids.txt --out .\results\kaks\target.cds.fa
```

脚本日志中的`Found`应与请求ID数量一致，`Missing`应为空。若没有全量CDS FASTA、只有基因组FASTA和GFF/GFF3，则本机原生Windows路线缺少简单、标准的gffread安装方式；请在Linux/WSL环境运行第9节的gffread命令生成全量CDS，再用与Windows路线一致的ID清单提取候选序列。不要使用不同版本的基因组和注释，也不要仅通过删改ID后缀来掩盖映射错误。

查看候选序列数和标题，并运行CDS长度检查：

```powershell
(Select-String -Path .\results\kaks\target.cds.fa -Pattern '^>').Count
Get-Content .\data\intermediate\kaks_ids.txt -TotalCount 10
py -3 .\scripts\validate_cds_lengths.py --cds .\results\kaks\target.cds.fa | Set-Content -Encoding UTF8 .\results\kaks\cds_validation.tsv
Get-Content .\results\kaks\cds_validation.tsv -TotalCount 10
```

`cds_validation.tsv`是供人工阅读的检查记录，不是Ka/Ks计算器输入。PowerShell中的`Set-Content -Encoding UTF8`在Windows PowerShell5.1中可能写入UTF-8 BOM；该检查表无需作为序列输入继续解析。按照第10节逐对复核CDS长度、阅读框、模糊碱基、遗传密码表、异构体及序列对应关系。

### 15.8 在TBtools-II中计算Ka/Ks

在TBtools-II中搜索并打开`Simple Ka/Ks Calculator (NG)`，按以下字段选择文件：`Input CDS File`选择`results\kaks\target.cds.fa`；`Input GenePair File`选择无表头、两列Tab分隔的`data\intermediate\kaks_pairs.tsv`；`Output Table File`设为`results\kaks\kaks.result.tsv`。输入框中应显示当前Windows项目路径；逐项核对后再运行。若界面提示ID缺失或输入格式错误，应先修复文件，不能将缺失配对解释为Ka/Ks为0。

计算完成后检查输出：

```powershell
Get-Item .\results\kaks\kaks.result.tsv | Select-Object FullName, Length
Get-Content .\results\kaks\kaks.result.tsv -TotalCount 10
Select-String -Path .\results\kaks\kaks.result.tsv -Pattern 'nan|inf' -CaseSensitive:$false
```

检查输出配对数与输入配对数是否相符，并逐项复核NaN、Inf、缺失行及异常Ka或Ks。解释原则和报告要求见第12节。

### 15.9 Windows结果计数、常见故障与提交检查

以下PowerShell代码按区分大小写的ID统计蛋白质FASTA唯一ID、位置文件ID和`.collinearity`中出现的唯一基因ID；它假定MCScanX配对行遵循标准`.collinearity`结构。运行后应将数量与脚本日志和原始文件人工抽查结果对照：

```powershell
$proteinIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
Get-Content .\data\intermediate\TK.clean.protein.fasta | ForEach-Object { if ($_ -match '^>(\S+)') { [void]$proteinIds.Add($Matches[1]) } }
$positionIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
Get-Content .\data\intermediate\TK.gff | ForEach-Object { $fields = $_.Split([char[]]@([char]9)); if ($fields.Length -ge 2) { [void]$positionIds.Add($fields[1]) } }
$collinearIds = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
Get-Content .\results\mcscanx\TK.collinearity | ForEach-Object {
    if ($_ -match '^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)') { [void]$collinearIds.Add($Matches[1]); [void]$collinearIds.Add($Matches[2]) }
}
$ratio = if ($positionIds.Count -gt 0) { [math]::Round(100 * $collinearIds.Count / $positionIds.Count, 2) } else { $null }
[pscustomobject]@{ ProteinUniqueIDs = $proteinIds.Count; PositionUniqueIDs = $positionIds.Count; CollinearUniqueGenes = $collinearIds.Count; CollinearPercentOfPositionIDs = $ratio }
```

若PowerShell提示`py`或`python`找不到命令，检查Python是否安装、`py -0p`是否列出解释器，并按实际可用启动器调整；不要把Conda Bash初始化命令粘贴进PowerShell。若DIAMOND提示运行库缺失，安装Microsoft Visual C++ Redistributable并重新验证。若找不到`.\tools\diamond\diamond.exe`，确认压缩包已解压且可执行文件位置准确，并使用`&`调用。若Windows版TBtools-II缺少`Quick Run MCScanX Wrapper`，或不能稳定生成`.collinearity`，将该阶段转到Linux/WSL。若需要从基因组提取CDS但无法安装gffread，也转到Linux/WSL完成该步骤。路径中包含空格时为每个参数加引号；输入FASTA、GFF/GFF3、蛋白质ID与CDS ID必须来自同一版本。

分析结束后可在PowerShell中检查Git状态，确认没有真实序列、未公开数据或中间结果被意外暂存：

```powershell
git status --short
```

Windows命令行步骤主要依赖PowerShell5.1及.NET文本流。若改用PowerShell7，仍可按本节执行，但需在更换文本编码参数前检查生成文件是否含BOM；DIAMOND、Python脚本及TBtools-II版本应记录在分析记录中。

## 参考资料

本操作文档中的脚本参数和输出定义以本仓库代码为准：[`prepare_mcscanx_inputs.py`](https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline/blob/main/scripts/prepare_mcscanx_inputs.py)、[`filter_collinearity_targets.py`](https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline/blob/main/scripts/filter_collinearity_targets.py)、[`extract_cds_by_ids.py`](https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline/blob/main/scripts/extract_cds_by_ids.py)和[`validate_cds_lengths.py`](https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline/blob/main/scripts/validate_cds_lengths.py)；TBtools-II输入栏名称和文件路径示例可对照[本仓库README](https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline/blob/main/README.md)。DIAMOND Windows二进制文件及运行库要求见[DIAMOND官方安装说明](https://github.com/bbuchfink/diamond_docs/blob/master/2%20Installation.MD)；MCScanX命令行支持环境见[MCScanX上游说明](https://github.com/wyp1125/MCScanX)；gffread平台及用法见[gffread项目说明](https://github.com/gpertea/gffread)和[Bioconda平台列表](https://anaconda.org/bioconda/gffread)；TBtools-II安装包见[官方发布页](https://github.com/CJ-Chen/TBtools-II/releases)。生物学及算法细节另见[MCScanX论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336/)和[DIAMOND官方文档](https://github.com/bbuchfink/diamond)。
