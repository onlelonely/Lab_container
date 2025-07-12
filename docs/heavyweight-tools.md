# 重量級工具容器指南

**版本**: v2.0  
**更新日期**: 2025年7月12日  

## 📋 概覽

本指南詳細說明 Lab Container 中重量級生物資訊工具的容器化使用方法，包括 GATK、ColabFold、AlphaFold 等大型分析軟體的安裝、配置和執行方式。

## 🏗️ 重量級工具特點

### 為什麼需要容器化？
- **資源需求大**: 記憶體需求 8GB+，儲存空間 50GB+
- **依賴複雜**: 特定 Python/Java/CUDA 版本依賴
- **環境衝突**: 不同工具間可能存在版本衝突
- **安裝困難**: 編譯時間長，配置複雜

### 容器化優勢
- **環境隔離**: 每個工具獨立運行環境
- **即時可用**: 預先建置的映像檔
- **資源控制**: 精確控制 CPU、記憶體使用
- **版本穩定**: 固定版本依賴關係

## 🧬 支援的重量級工具

### 基因體分析工具
- **GATK** - 基因變異檢測分析套件
- **FreeBayes** - 貝氏變異檢測器
- **Strelka2** - 體細胞變異檢測
- **Mutect2** - 癌症變異檢測

### 蛋白質結構預測
- **ColabFold** - 快速蛋白質摺疊預測
- **AlphaFold** - 深度學習結構預測
- **ChimeraX** - 分子視覺化工具

### 轉錄體分析
- **STAR** - RNA-seq 快速比對器
- **StringTie** - 轉錄體組裝和定量
- **Cufflinks** - 轉錄體分析套件

## 🚀 容器管理系統

### 容器啟動腳本

```bash
#!/bin/bash
# scripts/run-heavy-tool.sh

TOOL_NAME="$1"
TOOL_VERSION="$2"
shift 2

# 工具配置映射
declare -A TOOL_CONFIGS=(
    ["gatk"]="broadinstitute/gatk:${TOOL_VERSION:-4.4.0.0}"
    ["colabfold"]="colabfold/colabfold:${TOOL_VERSION:-1.5.5}"
    ["star"]="quay.io/biocontainers/star:${TOOL_VERSION:-2.7.10a}"
)

# 資源限制配置
declare -A RESOURCE_LIMITS=(
    ["gatk"]="--memory=16g --cpus=8"
    ["colabfold"]="--memory=32g --cpus=16 --gpus=all"
    ["star"]="--memory=32g --cpus=16"
)

# 檢查工具是否支援
if [[ ! ${TOOL_CONFIGS[$TOOL_NAME]} ]]; then
    echo "錯誤: 不支援的工具 '$TOOL_NAME'"
    echo "支援的工具: ${!TOOL_CONFIGS[@]}"
    exit 1
fi

# 準備容器環境
CONTAINER_IMAGE="${TOOL_CONFIGS[$TOOL_NAME]}"
RESOURCE_OPTS="${RESOURCE_LIMITS[$TOOL_NAME]}"
VOLUME_MOUNTS="-v $(pwd)/data:/data -v $(pwd)/results:/results"

echo "啟動工具: $TOOL_NAME"
echo "容器映像: $CONTAINER_IMAGE"

# 執行容器
docker run --rm -it \
    $RESOURCE_OPTS \
    $VOLUME_MOUNTS \
    --name "lab_${TOOL_NAME}_$(date +%s)" \
    $CONTAINER_IMAGE \
    "$@"
```

## 🧪 GATK - 基因變異檢測

### 環境準備

```bash
# 拉取 GATK 映像檔
docker pull broadinstitute/gatk:4.4.0.0

# 準備參考基因組
mkdir -p data/references
wget -O data/references/hg38.fasta \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
```

### 變異檢測流程

```bash
#!/bin/bash
# scripts/gatk-variant-calling.sh

SAMPLE_NAME="$1"
INPUT_BAM="data/aligned/${SAMPLE_NAME}_sorted.bam"
REFERENCE="data/references/hg38.fasta"
OUTPUT_VCF="results/variants/${SAMPLE_NAME}_variants.vcf"

echo "開始 GATK 變異檢測: $SAMPLE_NAME"

# 1. 標記重複序列
docker run --rm \
    -v $(pwd):/workspace \
    -w /workspace \
    --memory=16g --cpus=8 \
    broadinstitute/gatk:4.4.0.0 \
    gatk MarkDuplicates \
    -I $INPUT_BAM \
    -O data/aligned/${SAMPLE_NAME}_marked_duplicates.bam \
    -M data/aligned/${SAMPLE_NAME}_duplicate_metrics.txt

# 2. 變異檢測
docker run --rm \
    -v $(pwd):/workspace \
    -w /workspace \
    --memory=16g --cpus=8 \
    broadinstitute/gatk:4.4.0.0 \
    gatk HaplotypeCaller \
    -I data/aligned/${SAMPLE_NAME}_marked_duplicates.bam \
    -R $REFERENCE \
    -O $OUTPUT_VCF

echo "GATK 變異檢測完成: $SAMPLE_NAME"
```

## 🧬 ColabFold - 蛋白質摺疊預測

### 單一蛋白質預測

```bash
#!/bin/bash
# scripts/colabfold-prediction.sh

PROTEIN_FASTA="$1"
OUTPUT_DIR="results/colabfold/$(basename $PROTEIN_FASTA .fasta)"

echo "開始 ColabFold 蛋白質結構預測"

# 建立輸出目錄
mkdir -p $OUTPUT_DIR

# 執行結構預測
docker run --rm \
    -v $(pwd):/workspace \
    -w /workspace \
    --memory=32g --cpus=16 --gpus=all \
    colabfold/colabfold:1.5.5 \
    colabfold_batch \
    --num-models 5 \
    $PROTEIN_FASTA \
    $OUTPUT_DIR

echo "ColabFold 預測完成"
```

## 🌟 STAR - RNA-seq 比對器

### 基因組索引建立

```bash
#!/bin/bash
# scripts/star-index.sh

GENOME_FASTA="data/references/genome.fasta"
ANNOTATION_GTF="data/references/annotation.gtf"
INDEX_DIR="data/references/star_index"

echo "建立 STAR 基因組索引"

mkdir -p $INDEX_DIR

docker run --rm \
    -v $(pwd):/workspace \
    -w /workspace \
    --memory=64g --cpus=32 \
    quay.io/biocontainers/star:2.7.10a \
    STAR \
    --runMode genomeGenerate \
    --genomeDir $INDEX_DIR \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $ANNOTATION_GTF \
    --runThreads 32

echo "STAR 索引建立完成"
```

---

**相關文檔**:
- [輕量級工具使用](./lightweight-tools.md)
- [效能優化指南](./performance.md)
- [工作流程範例](./workflows.md)