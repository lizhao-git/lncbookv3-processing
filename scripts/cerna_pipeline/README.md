# cerna — 可复现的 ceRNA 分析流程（Python 3）

从遗留 ceRNA 脚本（Python 2、路径硬编码）重构而来的 ceRNA（competing
endogenous RNA）分析流程。整合三个 miRNA 靶点预测工具（miRanda、RNAhybrid、
TargetScan）的预测结果，取交集并经结合位点精炼得到预测的 ceRNA 相互作用，再结合
实验验证证据与疾病注释。

## 分析方案

```
三工具预测 (miRanda / RNAhybrid / TargetScan)
   ↓  归一化交互表  miRNA target score energy start end
   ↓  三工具交集：共同的 miRNA–lncRNA 对
   ↓  结合位点精炼：TargetScan 位点落在 miRanda/RNAhybrid 位点内，
   ↓                且 miRanda 位点上的 6-nt 滑窗落在 RNAhybrid 位点内
   ↓  three_overlap / three_overlap_unique  (预测的 ceRNA)
   ↓  实验验证：LncRNAWiki + HGNC → gene–miRNA
   ↓  最终表：transcript → gene (GeneID LncID miRNA Score Energy Start End)
   ↓  疾病注释：HMDD → lncRNA–disease
```

## 脚本

| 脚本 | 阶段 | 说明 |
|---|---|---|
| `common.py` | 公共库 | 归一化交互表加载 |
| `parse_tools.py` | 归一化 | 三个工具输出 → 统一交互表 |
| `predict_ceRNA.py` | 预测 | 三工具交集 + 结合位点精炼 |
| `experiment_annotate.py` | 下游 | 实验验证 + 最终表 + 疾病注释 |
| `pipeline.py` | 编排 | 读 JSON manifest，串行驱动上述阶段 |
| `prepare_sequences.py` | 上游·数据准备 | GTF+参考基因组→lncRNA.fa；mature.fa→hsa miRNA+种子 |
| `predict_targets.py` | 上游·预测 | 运行 miRanda/RNAhybrid/TargetScan 并归一化输出 |
| `pipeline_predict.py` | 上游·编排 | 读预测 manifest，串行驱动 prepare→predict |

## 依赖

- Python 3（仅标准库）
- 分析阶段（`pipeline.py`）：三个预测工具的输出文件作为**输入**
- 预测阶段（`pipeline_predict.py`）：`miranda` 与 `RNAhybrid` 二进制（已打包进
  Dockerfile）、参考基因组 hg38 FASTA、miRBase `mature.fa`

## 上游预测（从序列到预测文件）

`pipeline_predict.py` 从成熟 miRNA 与 GTF 中的 RNA 序列出发，运行三个预测工具，
生成分析阶段所需的三个输入文件：

```bash
python pipeline_predict.py --manifest cerna_predict.example.json \
  --data-root /path/to/data --output-dir results/cerna_predict
```

产物（`results/cerna_predict/predictions/` 下）：

| 文件 | 格式 |
|---|---|
| `miranda/miRanda_interaction_out_select_v2.txt` | `>miRNA target score energy start end` |
| `rnahybrid/RNAhybrid_out.txt` | `miRNA target start end energy pvalue` |
| `targetscan/TargetScan_combine_v2.txt` | `target miRNA start end` |

把这三项填进分析 manifest（`cerna.example.json`）的 `miranda`/`targetscan`/
`rnahybrid` 键即可衔接。预测 manifest 结构见
[`cerna_predict.example.json`](cerna_predict.example.json)：`genome_fa`、`gtf`、
`mature_mirna_fa` 按 `--data-root` 解析，另有工具与阈值项（`miranda_bin`、
`rnahybrid_bin`、`miranda_score`、`miranda_energy`、`rnahybrid_pvalue`、
`rnahybrid_model`）。

## 用法（直接调用）

```bash
python pipeline.py --manifest cerna.example.json \
  --data-root /path/to/data --output-dir results/cerna
```

`manifest` 结构见 [`cerna.example.json`](cerna.example.json)。

## 用法（Nextflow/Docker）

```bash
docker build --target cerna-predict -t lncbookv3-cerna-predict:latest .
nextflow run . --run_cerna true --cerna_manifest scripts/cerna_pipeline/cerna.example.json
```

单独运行上游预测工作流：

```bash
nextflow run ./subworkflows/local/cerna_predict/main.nf --manifest scripts/cerna_pipeline/cerna_predict.example.json
```

在顶层工作流中开启该分支：

```yaml
run_cerna: true
cerna_manifest: {class: File, path: path/to/cerna.example.json}
cerna_data_root: {class: Directory, path: path/to/data}
```

## 输入文件（manifest 相对路径）

| 键 | 说明 |
|---|---|
| `miranda` | `>miRNA target score energy start end` |
| `targetscan` | `target miRNA start end` |
| `rnahybrid` | `miRNA target start end energy pvalue`（目录或单文件）|
| `transcript_id_alias` | `gene_id ENST...;NR_...;XR_...` |
| `hgnc_symbol_transcript` | `symbol alias(|-sep) refseq` |
| `transcript_gene` | `transcript gene` |
| `lncrnawiki` | `lncRNA_symbol miRNA;miRNA;...` |
| `hmdd` | `... miRNA disease` |
