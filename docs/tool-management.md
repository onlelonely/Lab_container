# 工具管理系統指南

**文檔版本**: v2.0  
**核心組件**: Tool Management System  
**適用讀者**: 生物資訊研究者、進階使用者

## 🧬 工具管理系統概述

Lab Container 提供**混合工具管理架構**，將生物資訊工具分為**輕量級**（容器內整合）和**重量級**（專用容器）兩類，實現資源優化和使用便利的平衡。

### 設計理念
- **按需分配**: 根據工具特性選擇最適合的部署方式
- **資源優化**: 避免不必要的資源浪費
- **使用簡化**: 統一的命令介面管理所有工具
- **擴展性**: 支援新工具的輕鬆整合

## 🏗 工具分類架構

### 輕量級工具 (In-Container Tools)
**特徵**: 檔案小、啟動快、高頻使用
**部署方式**: 直接整合在主開發容器內

```bash
🔧 輕量級工具清單:
├── FastQC          # NGS 資料品質控制
├── SAMtools        # SAM/BAM 檔案處理
├── BCFtools        # VCF 檔案處理
├── BEDtools        # 基因組區間操作
├── BWA             # 序列比對工具
└── HTSlib          # 高通量測序資料處理庫
```

### 重量級工具 (Dedicated Container Tools)
**特徵**: 檔案大、資源需求高、專業化、低頻使用
**部署方式**: 獨立容器按需啟動

```bash
🚀 重量級工具清單:
├── GATK            # 基因組分析工具包
├── ColabFold       # 蛋白質結構預測
├── AutoDock Vina   # 分子對接模擬
├── BLAST+          # 序列相似性搜尋
└── 可擴展...        # 支援新增更多工具
```

## 🛠 工具管理器使用指南

### 基本命令格式
```bash
# 工具管理器位置
bash .devcontainer/scripts/manage/tool-manager.sh <操作> <工具名稱> [參數]
```

### 快速開始
```bash
# 顯示說明
bash .devcontainer/scripts/manage/tool-manager.sh help

# 檢查所有工具狀態
bash .devcontainer/scripts/manage/tool-manager.sh check all

# 管理工具容器
bash .devcontainer/scripts/manage/tool-manager.sh manage status
```

## 🔧 輕量級工具使用

### FastQC - NGS 品質控制
```bash
# 基本用法 (直接使用，無需容器)
fastqc input.fastq.gz

# 批次處理多個檔案
fastqc data/*.fastq.gz -o results/qc/

# 指定輸出目錄和格式
fastqc input.fastq.gz -o results/ --format=html

# 安靜模式執行
fastqc input.fastq.gz -q -o results/

# 檢查版本
fastqc --version
```

### SAMtools - SAM/BAM 處理
```bash
# 查看 BAM 檔案資訊
samtools view -H input.bam  # 查看標頭
samtools flagstat input.bam  # 統計資訊

# BAM 檔案操作
samtools sort input.bam -o sorted.bam  # 排序
samtools index sorted.bam  # 建立索引
samtools merge merged.bam file1.bam file2.bam  # 合併

# 格式轉換
samtools view -b input.sam > output.bam  # SAM 轉 BAM
samtools view -h input.bam > output.sam  # BAM 轉 SAM

# 過濾操作
samtools view -q 30 input.bam > high_quality.bam  # 品質過濾
samtools view -f 2 input.bam > paired.bam  # 配對reads

# 統計分析
samtools depth input.bam > depth.txt  # 深度統計
samtools coverage input.bam  # 覆蓋率統計
```

### BCFtools - VCF 處理
```bash
# VCF 檔案資訊
bcftools view -h variants.vcf  # 查看標頭
bcftools stats variants.vcf  # 統計資訊

# 檔案操作
bcftools sort variants.vcf -o sorted.vcf  # 排序
bcftools index sorted.vcf.gz  # 建立索引
bcftools merge file1.vcf file2.vcf -o merged.vcf  # 合併

# 過濾操作
bcftools filter -i 'QUAL>30' variants.vcf  # 品質過濾
bcftools view -m2 -M2 variants.vcf  # 雙等位基因位點

# 註釋操作
bcftools annotate -a annotations.bed variants.vcf
```

### BEDtools - 基因組區間操作
```bash
# 區間操作
bedtools intersect -a file1.bed -b file2.bed  # 交集
bedtools subtract -a file1.bed -b file2.bed   # 差集
bedtools merge -i sorted.bed  # 合併重疊區間

# 覆蓋率分析
bedtools coverage -a targets.bed -b alignments.bam

# 隨機取樣
bedtools sample -i input.bed -n 1000  # 隨機取樣1000個區間

# 格式轉換
bedtools bamtobed -i input.bam > output.bed
```

## 🚀 重量級工具使用

### GATK - 基因組分析工具包

#### 基本用法
```bash
# 檢查 GATK 可用性
bash .devcontainer/scripts/manage/tool-manager.sh check gatk

# 啟動 GATK 容器
bash .devcontainer/scripts/manage/tool-manager.sh manage start gatk

# 執行 GATK 命令
bash .devcontainer/scripts/manage/tool-manager.sh gatk --version
```

#### 常用 GATK 流程
```bash
# 1. 變異呼叫 (Variant Calling)
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I input.bam \
    -O variants.vcf \
    -R reference.fa

# 2. 變異品質分數重估 (VQSR)
bash .devcontainer/scripts/manage/tool-manager.sh gatk VariantRecalibrator \
    -V input.vcf \
    -O recal.table \
    --resource:dbsnp,known=true,training=false,truth=false dbsnp.vcf

# 3. 基因型呼叫
bash .devcontainer/scripts/manage/tool-manager.sh gatk GenotypeGVCFs \
    -V input.g.vcf \
    -O output.vcf \
    -R reference.fa

# 4. 變異過濾
bash .devcontainer/scripts/manage/tool-manager.sh gatk VariantFiltration \
    -V input.vcf \
    -O filtered.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "QD2"
```

### ColabFold - 蛋白質結構預測

#### 基本用法
```bash
# 檢查 ColabFold 可用性和 GPU 支援
bash .devcontainer/scripts/manage/tool-manager.sh check colabfold

# 單一蛋白質結構預測
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta \
    results/structure/

# 批次結構預測
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    proteins/ \
    results/batch_structure/ \
    --num-models 5
```

#### 進階選項
```bash
# 自訂預測參數
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta \
    results/ \
    --num-models 5 \
    --num-recycles 3 \
    --stop-at-score 90 \
    --use-gpu-relax

# MSA 模式選擇
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta \
    results/ \
    --msa-mode mmseqs2_uniref_env

# 模板使用
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta \
    results/ \
    --templates \
    --custom-template template.pdb
```

### AutoDock Vina - 分子對接

#### 基本用法
```bash
# 檢查 AutoDock Vina 可用性
bash .devcontainer/scripts/manage/tool-manager.sh check autodock

# 基本分子對接
bash .devcontainer/scripts/manage/tool-manager.sh autodock \
    --receptor receptor.pdbqt \
    --ligand ligand.pdbqt \
    --out docking_result.pdbqt \
    --log docking.log
```

#### 進階對接設置
```bash
# 指定對接範圍
bash .devcontainer/scripts/manage/tool-manager.sh autodock \
    --receptor receptor.pdbqt \
    --ligand ligand.pdbqt \
    --center_x 15.0 \
    --center_y 10.0 \
    --center_z 5.0 \
    --size_x 20 \
    --size_y 20 \
    --size_z 20 \
    --out result.pdbqt

# 多構象對接
bash .devcontainer/scripts/manage/tool-manager.sh autodock \
    --receptor receptor.pdbqt \
    --ligand ligand.pdbqt \
    --num_modes 10 \
    --energy_range 5 \
    --exhaustiveness 16 \
    --out multi_conf.pdbqt
```

## 🔄 容器管理操作

### 容器生命週期管理
```bash
# 啟動特定工具容器
bash .devcontainer/scripts/manage/tool-manager.sh manage start gatk
bash .devcontainer/scripts/manage/tool-manager.sh manage start colabfold
bash .devcontainer/scripts/manage/tool-manager.sh manage start autodock

# 停止工具容器
bash .devcontainer/scripts/manage/tool-manager.sh manage stop gatk
bash .devcontainer/scripts/manage/tool-manager.sh manage stop all

# 重啟容器
bash .devcontainer/scripts/manage/tool-manager.sh manage restart gatk

# 檢查容器狀態
bash .devcontainer/scripts/manage/tool-manager.sh manage status
```

### 資源管理
```bash
# 查看容器資源使用
docker stats

# 清理未使用的容器
bash .devcontainer/scripts/manage/tool-manager.sh manage clean

# 檢查磁碟使用
docker system df

# 清理所有未使用資源 (謹慎使用)
docker system prune -f
```

### 日誌管理
```bash
# 查看容器日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk

# 即時監控日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk -f

# 查看最近的日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk --tail 100
```

## 📊 工具狀態監控

### 健康狀態檢查
```bash
# 檢查所有工具可用性
bash .devcontainer/scripts/manage/tool-manager.sh check all

# 檢查特定工具
bash .devcontainer/scripts/manage/tool-manager.sh check gatk
bash .devcontainer/scripts/manage/tool-manager.sh check colabfold
bash .devcontainer/scripts/manage/tool-manager.sh check autodock

# 詳細狀態報告
bash .devcontainer/scripts/manage/devcontainer-status.sh
```

### 效能監控
```bash
# 容器資源使用情況
docker stats --no-stream

# 工具執行時間監控
time bash .devcontainer/scripts/manage/tool-manager.sh gatk <command>

# 磁碟空間監控
df -h
du -sh results/ data/ models/
```

## 🛠 工作流程整合

### NGS 資料分析完整流程
```bash
#!/bin/bash
# ngs_pipeline.sh - NGS 資料分析完整流程

# 設定變數
INPUT_DIR="data/fastq"
OUTPUT_DIR="results"
REFERENCE="data/reference/genome.fa"

echo "🧬 開始 NGS 資料分析流程"

# 1. 品質控制 (輕量級工具)
echo "📊 步驟 1: 品質控制"
mkdir -p "$OUTPUT_DIR/qc"
fastqc "$INPUT_DIR"/*.fastq.gz -o "$OUTPUT_DIR/qc"

# 2. 序列比對
echo "🎯 步驟 2: 序列比對"
mkdir -p "$OUTPUT_DIR/alignment"
for fastq in "$INPUT_DIR"/*.fastq.gz; do
    base=$(basename "$fastq" .fastq.gz)
    bwa mem "$REFERENCE" "$fastq" | \
    samtools sort -o "$OUTPUT_DIR/alignment/${base}.bam"
    samtools index "$OUTPUT_DIR/alignment/${base}.bam"
done

# 3. 變異呼叫 (重量級工具)
echo "🔬 步驟 3: 變異呼叫"
mkdir -p "$OUTPUT_DIR/variants"
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I "$OUTPUT_DIR/alignment/sample.bam" \
    -O "$OUTPUT_DIR/variants/raw_variants.vcf" \
    -R "$REFERENCE"

# 4. 變異過濾
echo "🔍 步驟 4: 變異過濾"
bash .devcontainer/scripts/manage/tool-manager.sh gatk VariantFiltration \
    -V "$OUTPUT_DIR/variants/raw_variants.vcf" \
    -O "$OUTPUT_DIR/variants/filtered_variants.vcf" \
    --filter-expression "QD < 2.0 || FS > 60.0" \
    --filter-name "basic_filter"

# 5. 結果統計
echo "📈 步驟 5: 結果統計"
bcftools stats "$OUTPUT_DIR/variants/filtered_variants.vcf" > \
    "$OUTPUT_DIR/variants/variant_stats.txt"

echo "✅ NGS 分析流程完成!"
```

### 蛋白質結構預測流程
```bash
#!/bin/bash
# structure_prediction.sh - 蛋白質結構預測流程

PROTEIN_FASTA="$1"
OUTPUT_DIR="$2"

if [ -z "$PROTEIN_FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "用法: $0 <protein.fasta> <output_directory>"
    exit 1
fi

echo "🧬 開始蛋白質結構預測"

# 1. 檢查輸入檔案
if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "❌ 找不到輸入檔案: $PROTEIN_FASTA"
    exit 1
fi

# 2. 建立輸出目錄
mkdir -p "$OUTPUT_DIR"

# 3. 檢查 GPU 可用性
echo "🎮 檢查 GPU 狀態..."
if nvidia-smi &>/dev/null; then
    echo "✅ NVIDIA GPU 可用"
    GPU_FLAG="--use-gpu-relax"
else
    echo "⚠️  未檢測到 GPU，將使用 CPU 模式"
    GPU_FLAG=""
fi

# 4. 執行結構預測
echo "🔬 執行 ColabFold 結構預測..."
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    "$PROTEIN_FASTA" \
    "$OUTPUT_DIR" \
    --num-models 5 \
    --num-recycles 3 \
    $GPU_FLAG

# 5. 結果整理
echo "📊 整理預測結果..."
find "$OUTPUT_DIR" -name "*.pdb" | while read pdb_file; do
    echo "PDB 檔案: $pdb_file"
    echo "檔案大小: $(du -h "$pdb_file" | cut -f1)"
done

echo "✅ 蛋白質結構預測完成!"
echo "結果位置: $OUTPUT_DIR"
```

## ⚙️ 自訂工具整合

### 新增輕量級工具
```bash
# 1. 更新 Profile 配置
# 編輯 .devcontainer/profiles/bioinformatics.env
BIOINFORMATICS_TOOLS="fastqc samtools bcftools bedtools newtool"

# 2. 更新 Conda 環境檔案
# 編輯 .devcontainer/configs/conda/environment-bioinformatics.yml
dependencies:
  - bioconda::newtool

# 3. 重新生成環境
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics
```

### 新增重量級工具
```bash
# 1. 更新 docker-compose.tools.yml
# 新增服務定義
  newtool:
    image: newtool/newtool:latest
    volumes:
      - ./data:/workspace/data
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["custom", "newtool"]

# 2. 更新工具管理腳本
# 編輯 .devcontainer/scripts/manage/tool-manager.sh
# 新增工具包裝函數
run_newtool() {
    echo "🔧 執行 NewTool: $@"
    docker-compose --profile newtool run --rm newtool "$@"
}

# 3. 測試新工具
bash .devcontainer/scripts/manage/tool-manager.sh newtool --version
```

## 🚨 故障排除

### 輕量級工具問題
```bash
# 問題: 工具找不到
# 解決方案: 檢查是否正確安裝
which fastqc
conda list | grep fastqc

# 重新安裝
conda install -c bioconda fastqc

# 問題: 權限錯誤
# 解決方案: 檢查檔案權限
ls -la input.fastq.gz
chmod +r input.fastq.gz
```

### 重量級工具問題
```bash
# 問題: 容器啟動失敗
# 解決方案: 檢查 Docker 狀態
docker ps -a
docker logs <container_name>

# 重新拉取映像
docker pull broadinstitute/gatk:latest

# 問題: GPU 不可用
# 解決方案: 檢查 NVIDIA 驅動
nvidia-smi
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi

# 問題: 記憶體不足
# 解決方案: 調整容器資源限制
# 編輯 docker-compose.tools.yml
deploy:
  resources:
    limits:
      memory: 16G
```

### 網路連線問題
```bash
# 問題: 無法下載工具映像
# 解決方案: 檢查網路連線
ping hub.docker.com

# 使用國內映像源
docker pull registry.cn-hangzhou.aliyuncs.com/tool/gatk:latest
```

## 🎯 最佳實踐

### 1. 資源管理
- **按需啟動**: 只啟動需要的工具容器
- **及時清理**: 定期清理不需要的容器和映像
- **監控資源**: 定期檢查磁碟和記憶體使用

### 2. 資料管理
- **路徑一致**: 使用統一的資料目錄結構
- **備份重要**: 及時備份重要的分析結果
- **版本記錄**: 記錄使用的工具版本

### 3. 工作流程
- **腳本化**: 將常用分析流程寫成腳本
- **模組化**: 將複雜流程拆分為可重用的模組
- **文檔化**: 詳細記錄分析步驟和參數

---

**透過統一的工具管理系統，讓生物資訊分析變得更加簡單高效**  
*輕量級與重量級工具的完美結合，滿足各種分析需求！*