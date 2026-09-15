# methylation — 可复现的甲基化分析流程（Python 3）

从 `meth/` 遗留脚本（Python 2、路径硬编码、列索引/坐标体系不一致）重构而来，
并补写了缺失的阶段 B 差异检验。

> 与遗留脚本不同，本流程**不再区分 lncRNA/PCG**，而是基于通用 GTF 推导的基因列表
> （manifest 的 `gene_list`，可选；缺省时分析差异检验产出的全部基因）进行统计。

这些脚本被打包进 Docker 镜像（`COPY scripts/` → `/opt/lncbookv3/scripts/`），
由 Nextflow 分支 `subworkflows/local/methylation/main.nf` 调用。

## 脚本

| 脚本 | 阶段 | 说明 |
|---|---|---|
| `common.py` | 公共库 | 基因信息加载、区间交叠求均值、BH-FDR、曼-惠特尼 U（纯标准库）|
| `preprocess.py` | A1–A3 | 三种格式统一为 4 列 bed + 可选 hg19→hg38 liftover |
| `extract.py` | A4–A5 | 基因体/启动子区域平均甲基化 → 基因×样本矩阵 |
| `difftest.py` | B | Wilcoxon+FDR / 一致性+倍数变化 差异检验 |
| `aggregate.py` | C | 显著性矩阵、≥5&0 计数筛选、标签、数据库导入矩阵 |
| `beta_convert.py` | D | 百分比→beta 值（修复缩进 bug）|
| `pipeline.py` | 编排 | 读 JSON manifest，串行驱动上述阶段 |

## 依赖

- Python 3（仅标准库）
- `preprocess.py` 的外部工具（按格式）：
  - `bigWigToWig` + `wig2bed`（BEDOPS）— bigwig
  - `bigBedToBed` — bigbed
  - `liftOver` — 需要坐标转换时（UCSC kent tools）
  - bismark / bed 格式无需外部工具

## 用法（直接调用）

```bash
python pipeline.py --manifest methylation.example.json \
  --data-root /path/to/data --output-dir results/methylation
```

`manifest` 结构见 [`methylation.example.json`](methylation.example.json)。

## 用法（Nextflow/Docker）

```bash
docker build --target kent -t lncbookv3-kent:latest .
nextflow run . --run_methylation true --methylation_manifest scripts/methylation/methylation.example.json
```

在顶层工作流中开启该分支：

```yaml
run_methylation: true
methylation_manifest: {class: File, path: path/to/methylation.example.json}
methylation_data_root: {class: Directory, path: path/to/data}
```

## 修复的遗留问题

| 遗留问题 | 重构后的处理 |
|---|---|
| 路径硬编码 `/disk1/...` | 全部通过 manifest / CLI 参数传入 |
| Python 2 | 全部 Python 3（仅标准库）|
| 值列不一致 | 统一输出 4 列 bed，值恒为第 4 列 |
| 0/1 索引混用 | 统一 0-based 半开区间 `[start, end)` |
| step5 按 1 万基因分块 | 单遍处理，每样本每染色体只加载一次 |
| `beta_value.py` 缩进 bug | 已修复 |
| 阶段 B 脚本缺失 | 补写 `difftest.py` |
