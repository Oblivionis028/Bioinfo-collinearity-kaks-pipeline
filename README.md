<h1 align="center">共线性分析与Ka/Ks分析流程</h1>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.x-3776AB?style=flat-square&logo=python&logoColor=white" alt="Python 3" />
  <img src="https://img.shields.io/badge/Shell-Bash-4EAA25?style=flat-square&logo=gnubash&logoColor=white" alt="Bash" />
  <img src="https://img.shields.io/badge/MCScanX-Collinearity-2E8B57?style=flat-square" alt="MCScanX" />
  <img src="https://img.shields.io/badge/Ka%2FKs-Analysis-DC143C?style=flat-square" alt="Ka/Ks" />
  <img src="https://img.shields.io/badge/Data-Not%20Included-lightgrey?style=flat-square" alt="Data not included" />
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" alt="MIT License" />
</p>

本仓库提供可复用的种内或种间共线性分析与Ka/Ks分析模板，适用于从基因组注释、蛋白质序列和CDS序列出发，准备MCScanX输入，筛选目标基因相关的共线性基因对，提取候选CDS并开展基础质量检查。流程说明、脚本和示例配置均不包含真实基因组、注释、序列或分析结果。

## 项目功能

- 清理蛋白质FASTA标题，并从GFF/GFF3注释中提取坐标以生成MCScanX四列GFF文件。
- 配合DIAMOND开展蛋白质全对全比对，并通过MCScanX或TBtools-II识别共线性区块。
- 从`.collinearity`文件筛选与目标基因相关的共线性基因对，并从全量CDS FASTA提取对应CDS。
- 检查CDS长度、密码子相关特征和非法字符，随后配合TBtools-II计算Ka、Ks及Ka/Ks。

## 流程边界

本仓库不是一键式全自动流水线。Python脚本负责数据准备、筛选和基础质量检查；DIAMOND、MCScanX和Ka/Ks计算由外部软件完成。流程不会自动下载参考数据、安装第三方软件、判断GFF应使用哪个ID字段、自动建立基因与转录本及蛋白质ID映射，也不会替代对同源关系、基因模型、序列比对和异常结果的人工检查。

共线性关系可为基因组区域及基因复制历史提供证据，但不等同于同源关系的最终判定，也不能单独证明一对基因是直系同源基因。Ka/Ks是针对一对同源编码序列、在特定序列比对和替代模型下估算的比值；单个基因对的Ka/Ks大于1本身不足以证明适应性进化。详细的生物学注意事项、输入要求、命令和结果判读见[操作文档](docs/workflow.zh-CN.md)。

## 环境与依赖

推荐使用Linux/WSL完成完整命令行流程；Windows原生环境也可运行本仓库Python脚本、DIAMOND Windows程序及TBtools-II图形界面，但MCScanX封装是否可用取决于TBtools-II Windows版本。若需从基因组FASTA和GFF/GFF3生成CDS，推荐在Linux/WSL中使用[gffread](https://github.com/gpertea/gffread)。Python脚本仅使用标准库，不需要额外安装Python包。必需软件为Python 3、[DIAMOND](https://github.com/bbuchfink/diamond)和[TBtools-II](https://github.com/CJ-Chen/TBtools-II)；完整安装说明和两套逐步操作路线见[操作文档](docs/workflow.zh-CN.md)。

以下是Linux/WSL中的Bash检查命令；Windows原生PowerShell检查方式见[操作文档第15节](docs/workflow.zh-CN.md#十五windows原生环境powershell操作路线)：

```bash
python3 --version
diamond version
gffread --version    # 仅在需要时检查
```

TBtools-II为图形界面软件，需单独启动并确认能够正常运行。

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
│   ├── raw/
│   └── intermediate/
├── results/
│   ├── mcscanx/
│   └── kaks/
├── docs/
│   └── workflow.zh-CN.md
├── .gitignore
├── LICENSE
└── README.md
```

`data/`和`results/`中的实际数据与分析结果默认由`.gitignore`排除。提交前仍须检查Git暂存区，避免意外上传真实序列、未公开数据或样本信息。

## 快速开始

```bash
git clone https://github.com/Oblivionis028/Bioinfo-collinearity-kaks-pipeline.git
cd Bioinfo-collinearity-kaks-pipeline
```

将输入文件放入`data/raw/`，并先根据实际GFF/GFF3确认用于匹配蛋白质FASTA的feature类型和属性字段。完整步骤及参数解释见[操作文档](docs/workflow.zh-CN.md)。

## 主要结果文件

| 文件 | 用途 |
|---|---|
| `TK.clean.protein.fasta` | 标题清理后的蛋白质FASTA |
| `TK.gff` | MCScanX四列坐标文件 |
| `TK.blast` | DIAMOND蛋白质比对结果 |
| `TK.collinearity` | MCScanX共线性区块结果 |
| `target_collinearity_pairs.tsv` | 目标相关共线性基因对及坐标 |
| `target_collinearity_summary.tsv` | 目标基因命中情况汇总 |
| `target.cds.fa` | 候选基因对CDS序列 |
| `cds_validation.tsv` | CDS基础质量检查结果 |
| `kaks.result.tsv` | Ka、Ks及Ka/Ks计算结果 |

## 数据与许可证

仓库默认忽略FASTA、GFF/GTF、GenBank文件、比对数据库与结果、共线性结果以及常见表格、图片、日志和报告文件。上传前请运行`git status`并检查待提交内容。本项目采用[MIT License](LICENSE)。
