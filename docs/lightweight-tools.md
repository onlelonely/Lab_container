# 輕量級工具使用指南

**版本**: v2.0  
**更新日期**: 2025年7月12日  

## 📋 概覽

本指南詳細說明 Lab Container 中整合的輕量級生物資訊工具使用方法，包括品質控制、序列處理、比對分析等常用工具。

## 🧬 支援的輕量級工具

### 序列品質控制
- **FastQC** - 高通量測序資料品質評估
- **MultiQC** - 多樣本品質報告彙整
- **Trimmomatic** - 序列修剪和過濾
- **Cutadapt** - 接頭序列移除

### 序列比對工具
- **BWA** - Burrows-Wheeler 比對器
- **Bowtie2** - 快速序列比對工具
- **Minimap2** - 長讀長序列比對
- **HISAT2** - RNA-seq 比對工具

### 檔案處理工具
- **SAMtools** - SAM/BAM 檔案處理套件
- **BEDtools** - 基因組區間運算工具
- **BCFtools** - VCF/BCF 檔案處理
- **Picard** - 高通量測序資料處理

## 🚀 快速開始

### 環境啟動

```bash
# 啟動包含輕量級工具的環境
./scripts/start-lightweight-env.sh

# 檢查工具安裝狀態
./scripts/check-tools.sh lightweight
```

## 🔍 品質控制工具

### FastQC - 序列品質評估

```bash
# 單一檔案分析
fastqc data/raw/sample1_R1.fastq.gz \
    --outdir data/qc/ \
    --format fastq \
    --threads 4

# 批次分析
fastqc data/raw/*.fastq.gz \
    --outdir data/qc/ \
    --threads 8
```

### MultiQC - 多樣本報告彙整

```bash
# 彙整所有 FastQC 結果
multiqc data/qc/ \
    --outdir data/qc/multiqc/ \
    --title "Lab_Container_QC_Report"
```

## 🧭 序列比對工具

### BWA - 序列比對

```bash
# 建立 BWA 索引
bwa index data/genomes/reference.fasta

# 執行比對
bwa mem -t 8 \
    data/genomes/reference.fasta \
    data/trimmed/sample1_R1.fastq.gz \
    data/trimmed/sample1_R2.fastq.gz \
    > data/aligned/sample1.sam
```

### SAMtools - SAM/BAM 處理

```bash
# SAM 轉 BAM 並排序
samtools view -bS data/aligned/sample1.sam | \
    samtools sort -o data/aligned/sample1_sorted.bam

# 建立索引
samtools index data/aligned/sample1_sorted.bam

# 檢視比對統計
samtools flagstat data/aligned/sample1_sorted.bam
```

## 📊 批次處理工作流程

### RNA-seq 分析流程

```bash
#!/bin/bash
# scripts/rnaseq_pipeline.sh

SAMPLES_DIR="data/raw"
OUTPUT_DIR="data/processed"

# 建立輸出目錄
mkdir -p ${OUTPUT_DIR}/{qc,trimmed,aligned}

# 處理每個樣本
for sample in ${SAMPLES_DIR}/*_R1.fastq.gz; do
    base=$(basename $sample _R1.fastq.gz)
    echo "處理樣本: $base"
    
    # 品質控制
    fastqc ${SAMPLES_DIR}/${base}_R*.fastq.gz \
        --outdir ${OUTPUT_DIR}/qc \
        --threads 4
    
    # 序列比對
    bwa mem -t 8 data/genomes/reference.fasta \
        ${SAMPLES_DIR}/${base}_R1.fastq.gz \
        ${SAMPLES_DIR}/${base}_R2.fastq.gz \
        > ${OUTPUT_DIR}/aligned/${base}.sam
    
    # 轉換和排序
    samtools view -bS ${OUTPUT_DIR}/aligned/${base}.sam | \
        samtools sort -o ${OUTPUT_DIR}/aligned/${base}_sorted.bam
    
    echo "樣本 $base 處理完成"
done

echo "RNA-seq 流程完成"
```

---

**相關文檔**:
- [重量級工具容器](./heavyweight-tools.md)
- [工作流程範例](./workflows.md)
- [套件管理](./package-management.md)