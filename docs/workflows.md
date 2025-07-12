# 工作流程範例指南

**文檔版本**: v2.0  
**適用領域**: NGS分析、蛋白質結構預測、機器學習、資料科學  
**目標讀者**: 研究者、分析師、學習者

## 🧬 NGS 資料分析工作流程

### 基本 NGS 分析流程

#### 完整基因組重測序分析
```bash
#!/bin/bash
# wgs_analysis.sh - 全基因組重測序分析流程
# 用法: bash wgs_analysis.sh <sample_name> <fastq_dir> <reference_genome>

set -e  # 遇到錯誤即停止

SAMPLE_NAME="$1"
FASTQ_DIR="$2"
REFERENCE="$3"
OUTPUT_DIR="results/${SAMPLE_NAME}"

# 檢查參數
if [ $# -ne 3 ]; then
    echo "用法: $0 <sample_name> <fastq_dir> <reference_genome>"
    echo "範例: $0 sample01 data/fastq data/reference/hg38.fa"
    exit 1
fi

echo "🧬 開始全基因組重測序分析: $SAMPLE_NAME"
echo "================================================"

# 建立輸出目錄
mkdir -p "$OUTPUT_DIR"/{qc,alignment,variants,reports}

# ============================================================================
# 步驟 1: 原始資料品質控制
# ============================================================================
echo "📊 步驟 1/6: 原始資料品質控制"
echo "輸入目錄: $FASTQ_DIR"

# FastQC 品質檢查
fastqc "$FASTQ_DIR"/${SAMPLE_NAME}*.fastq.gz -o "$OUTPUT_DIR/qc/" -t 4

# MultiQC 整合報告 (如果可用)
if command -v multiqc &> /dev/null; then
    multiqc "$OUTPUT_DIR/qc/" -o "$OUTPUT_DIR/reports/" -n "${SAMPLE_NAME}_fastqc_report"
fi

echo "✅ 品質控制完成"

# ============================================================================  
# 步驟 2: 讀段比對
# ============================================================================
echo "🎯 步驟 2/6: 序列比對"

# 檢查參考基因組索引
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "建立 BWA 索引..."
    bwa index "$REFERENCE"
fi

# BWA 比對 (假設為配對末端測序)
READ1="$FASTQ_DIR/${SAMPLE_NAME}_R1.fastq.gz"
READ2="$FASTQ_DIR/${SAMPLE_NAME}_R2.fastq.gz"

if [ -f "$READ1" ] && [ -f "$READ2" ]; then
    echo "配對末端比對: $READ1, $READ2"
    bwa mem -t 8 -M "$REFERENCE" "$READ1" "$READ2" | \
    samtools view -Sb - | \
    samtools sort -@ 4 -o "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam"
else
    echo "❌ 找不到配對的 FASTQ 檔案"
    exit 1
fi

# 建立 BAM 索引
samtools index "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam"

echo "✅ 序列比對完成"

# ============================================================================
# 步驟 3: 比對後處理
# ============================================================================
echo "🔧 步驟 3/6: 比對後處理"

# 標記重複序列
samtools markdup \
    "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam" \
    "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam"

samtools index "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam"

# 比對統計
samtools flagstat "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam" > \
    "$OUTPUT_DIR/reports/${SAMPLE_NAME}.alignment_stats.txt"

samtools coverage "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam" > \
    "$OUTPUT_DIR/reports/${SAMPLE_NAME}.coverage_stats.txt"

echo "✅ 比對後處理完成"

# ============================================================================
# 步驟 4: 變異檢測 (使用 GATK)
# ============================================================================
echo "🔬 步驟 4/6: 變異檢測"

# 確保 GATK 容器可用
bash .devcontainer/scripts/manage/tool-manager.sh manage start gatk

# GATK HaplotypeCaller 變異呼叫
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam" \
    -O "$OUTPUT_DIR/variants/${SAMPLE_NAME}.raw.vcf.gz" \
    -R "$REFERENCE" \
    --emit-ref-confidence GVCF

echo "✅ 變異檢測完成"

# ============================================================================
# 步驟 5: 變異過濾
# ============================================================================
echo "🔍 步驟 5/6: 變異過濾"

# 基本品質過濾
bash .devcontainer/scripts/manage/tool-manager.sh gatk VariantFiltration \
    -V "$OUTPUT_DIR/variants/${SAMPLE_NAME}.raw.vcf.gz" \
    -O "$OUTPUT_DIR/variants/${SAMPLE_NAME}.filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "basic_filter"

# 進一步過濾 (移除低品質變異)
bcftools filter -i 'FILTER="PASS"' \
    "$OUTPUT_DIR/variants/${SAMPLE_NAME}.filtered.vcf.gz" \
    -o "$OUTPUT_DIR/variants/${SAMPLE_NAME}.final.vcf.gz"

echo "✅ 變異過濾完成"

# ============================================================================
# 步驟 6: 結果統計與報告
# ============================================================================
echo "📈 步驟 6/6: 生成分析報告"

# VCF 統計
bcftools stats "$OUTPUT_DIR/variants/${SAMPLE_NAME}.final.vcf.gz" > \
    "$OUTPUT_DIR/reports/${SAMPLE_NAME}.variant_stats.txt"

# 變異類型統計
bcftools query -f '%TYPE\n' "$OUTPUT_DIR/variants/${SAMPLE_NAME}.final.vcf.gz" | \
    sort | uniq -c > "$OUTPUT_DIR/reports/${SAMPLE_NAME}.variant_types.txt"

# 生成最終報告
cat > "$OUTPUT_DIR/reports/${SAMPLE_NAME}_final_report.txt" << EOF
========================================
NGS 分析報告: $SAMPLE_NAME
========================================
分析日期: $(date)
參考基因組: $REFERENCE

檔案位置:
- 比對檔案: $OUTPUT_DIR/alignment/${SAMPLE_NAME}.markdup.bam
- 變異檔案: $OUTPUT_DIR/variants/${SAMPLE_NAME}.final.vcf.gz
- 品質報告: $OUTPUT_DIR/qc/
- 統計報告: $OUTPUT_DIR/reports/

主要統計:
$(cat "$OUTPUT_DIR/reports/${SAMPLE_NAME}.alignment_stats.txt")

變異統計:
$(head -20 "$OUTPUT_DIR/reports/${SAMPLE_NAME}.variant_stats.txt")
========================================
EOF

echo "✅ NGS 分析流程完成!"
echo "📂 結果位置: $OUTPUT_DIR"
echo "📋 檢視報告: $OUTPUT_DIR/reports/${SAMPLE_NAME}_final_report.txt"
```

### RNA-seq 表達分析流程
```bash
#!/bin/bash
# rnaseq_analysis.sh - RNA-seq 表達分析流程

SAMPLE_NAME="$1"
FASTQ_DIR="$2"
REFERENCE_GENOME="$3"
ANNOTATION_GTF="$4"
OUTPUT_DIR="results/rnaseq/${SAMPLE_NAME}"

echo "🧬 開始 RNA-seq 表達分析: $SAMPLE_NAME"

mkdir -p "$OUTPUT_DIR"/{qc,alignment,counts,reports}

# 1. 品質控制
echo "📊 步驟 1: 品質控制"
fastqc "$FASTQ_DIR"/${SAMPLE_NAME}*.fastq.gz -o "$OUTPUT_DIR/qc/"

# 2. 序列比對 (使用 STAR 或 HISAT2)
echo "🎯 步驟 2: RNA-seq 比對"
if command -v hisat2 &> /dev/null; then
    # 使用 HISAT2
    hisat2 -x "$REFERENCE_GENOME" \
           -1 "$FASTQ_DIR/${SAMPLE_NAME}_R1.fastq.gz" \
           -2 "$FASTQ_DIR/${SAMPLE_NAME}_R2.fastq.gz" \
           -S "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sam" \
           --threads 8
    
    # 轉換為 BAM 並排序
    samtools view -bS "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sam" | \
    samtools sort -o "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam"
    
    samtools index "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam"
    rm "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sam"  # 節省空間
fi

# 3. 基因計數
echo "🔢 步驟 3: 基因表達計數"
if command -v featureCounts &> /dev/null; then
    featureCounts -a "$ANNOTATION_GTF" \
                  -o "$OUTPUT_DIR/counts/${SAMPLE_NAME}.counts.txt" \
                  -T 8 \
                  -p \
                  "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam"
fi

# 4. 表達統計
echo "📊 步驟 4: 表達統計"
samtools flagstat "$OUTPUT_DIR/alignment/${SAMPLE_NAME}.sorted.bam" > \
    "$OUTPUT_DIR/reports/${SAMPLE_NAME}.alignment_stats.txt"

echo "✅ RNA-seq 分析完成!"
```

## 🧪 蛋白質結構預測工作流程

### 單一蛋白質結構預測
```bash
#!/bin/bash
# protein_structure_prediction.sh - 蛋白質結構預測流程

PROTEIN_FASTA="$1"
PROTEIN_NAME="$2"
OUTPUT_DIR="results/structure/${PROTEIN_NAME}"

if [ $# -ne 2 ]; then
    echo "用法: $0 <protein.fasta> <protein_name>"
    exit 1
fi

echo "🧬 開始蛋白質結構預測: $PROTEIN_NAME"
echo "============================================"

mkdir -p "$OUTPUT_DIR"/{models,analysis,reports}

# ============================================================================
# 步驟 1: 序列驗證
# ============================================================================
echo "🔍 步驟 1/5: 序列驗證"

# 檢查 FASTA 格式
if ! grep -q "^>" "$PROTEIN_FASTA"; then
    echo "❌ 無效的 FASTA 格式"
    exit 1
fi

# 序列統計
SEQUENCE_LENGTH=$(grep -v "^>" "$PROTEIN_FASTA" | tr -d '\n' | wc -c)
echo "序列長度: $SEQUENCE_LENGTH 個氨基酸"

if [ "$SEQUENCE_LENGTH" -gt 2000 ]; then
    echo "⚠️  警告: 序列較長 (>2000 AA)，預測時間可能很長"
fi

echo "✅ 序列驗證完成"

# ============================================================================
# 步驟 2: GPU 環境檢查
# ============================================================================
echo "🎮 步驟 2/5: GPU 環境檢查"

if nvidia-smi &>/dev/null; then
    echo "✅ NVIDIA GPU 可用"
    GPU_MEMORY=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    echo "GPU 記憶體: ${GPU_MEMORY} MB"
    
    if [ "$GPU_MEMORY" -lt 8000 ]; then
        echo "⚠️  警告: GPU 記憶體可能不足，建議使用 CPU 模式"
        USE_GPU=false
    else
        USE_GPU=true
    fi
else
    echo "ℹ️  未檢測到 GPU，將使用 CPU 模式"
    USE_GPU=false
fi

# ============================================================================
# 步驟 3: ColabFold 結構預測
# ============================================================================
echo "🔬 步驟 3/5: ColabFold 結構預測"

# 確保 ColabFold 容器可用
bash .devcontainer/scripts/manage/tool-manager.sh check colabfold

# 設定預測參數
if [ "$USE_GPU" = true ]; then
    COLABFOLD_ARGS="--use-gpu-relax --num-models 5 --num-recycles 3"
else
    COLABFOLD_ARGS="--num-models 3 --num-recycles 1"
fi

# 執行結構預測
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    "$PROTEIN_FASTA" \
    "$OUTPUT_DIR/models/" \
    $COLABFOLD_ARGS

echo "✅ 結構預測完成"

# ============================================================================
# 步驟 4: 結果分析
# ============================================================================
echo "📊 步驟 4/5: 結果分析"

# 找到最佳模型 (通常是 rank_001)
BEST_MODEL=$(find "$OUTPUT_DIR/models/" -name "*rank_001*.pdb" | head -1)

if [ -n "$BEST_MODEL" ]; then
    echo "最佳模型: $BEST_MODEL"
    
    # 複製最佳模型到分析目錄
    cp "$BEST_MODEL" "$OUTPUT_DIR/analysis/${PROTEIN_NAME}_best_model.pdb"
    
    # 提取信心分數 (如果有 JSON 檔案)
    JSON_FILE="${BEST_MODEL%.pdb}.json"
    if [ -f "$JSON_FILE" ]; then
        python3 << EOF
import json
import numpy as np

with open('$JSON_FILE', 'r') as f:
    data = json.load(f)

if 'confidenceScore' in data:
    scores = data['confidenceScore']
    avg_score = np.mean(scores)
    print(f"平均信心分數: {avg_score:.2f}")
    
    high_conf = sum(1 for s in scores if s > 90)
    total = len(scores)
    print(f"高信心區域: {high_conf}/{total} ({high_conf/total*100:.1f}%)")
EOF
    fi
else
    echo "❌ 未找到預測模型"
fi

# 統計所有模型
MODEL_COUNT=$(find "$OUTPUT_DIR/models/" -name "*.pdb" | wc -l)
echo "生成模型數量: $MODEL_COUNT"

echo "✅ 結果分析完成"

# ============================================================================
# 步驟 5: 生成報告
# ============================================================================
echo "📋 步驟 5/5: 生成分析報告"

cat > "$OUTPUT_DIR/reports/${PROTEIN_NAME}_structure_report.txt" << EOF
========================================
蛋白質結構預測報告
========================================
蛋白質名稱: $PROTEIN_NAME
輸入檔案: $PROTEIN_FASTA
序列長度: $SEQUENCE_LENGTH 個氨基酸
預測日期: $(date)
使用 GPU: $USE_GPU

預測結果:
- 模型數量: $MODEL_COUNT
- 最佳模型: $(basename "$BEST_MODEL" 2>/dev/null || echo "無")
- 結果目錄: $OUTPUT_DIR/models/

建議後續分析:
1. 使用 PyMOL 或 ChimeraX 視覺化結構
2. 進行結構品質評估 (Ramachandran plot)
3. 功能域分析和活性位點預測
4. 與已知結構進行比較 (如有需要)

檔案說明:
- *.pdb: 蛋白質結構檔案
- *.json: 預測信心分數和元資料
- *_coverage.png: MSA 覆蓋率圖表 (如有)
- *_PAE.png: 預測校準誤差圖表 (如有)
========================================
EOF

echo "✅ 蛋白質結構預測流程完成!"
echo "📂 結果位置: $OUTPUT_DIR"
echo "📋 查看報告: $OUTPUT_DIR/reports/${PROTEIN_NAME}_structure_report.txt"

# 快速檢視最佳模型 (如果安裝了 PyMOL)
if command -v pymol &> /dev/null && [ -n "$BEST_MODEL" ]; then
    echo "🔬 是否要用 PyMOL 開啟結構? (y/N)"
    read -r response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        pymol "$BEST_MODEL" &
    fi
fi
```

### 批次蛋白質結構預測
```bash
#!/bin/bash
# batch_structure_prediction.sh - 批次蛋白質結構預測

FASTA_DIR="$1"
OUTPUT_BASE="results/batch_structure"

if [ ! -d "$FASTA_DIR" ]; then
    echo "用法: $0 <fasta_directory>"
    exit 1
fi

echo "🧬 開始批次蛋白質結構預測"
echo "輸入目錄: $FASTA_DIR"

mkdir -p "$OUTPUT_BASE"

# 找到所有 FASTA 檔案
FASTA_FILES=($(find "$FASTA_DIR" -name "*.fasta" -o -name "*.fa"))

if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "❌ 在 $FASTA_DIR 中未找到 FASTA 檔案"
    exit 1
fi

echo "找到 ${#FASTA_FILES[@]} 個 FASTA 檔案"

# 批次處理
for fasta_file in "${FASTA_FILES[@]}"; do
    protein_name=$(basename "$fasta_file" | sed 's/\.\(fasta\|fa\)$//')
    echo "處理: $protein_name"
    
    # 呼叫單一蛋白質預測腳本
    bash protein_structure_prediction.sh "$fasta_file" "$protein_name"
    
    echo "完成: $protein_name"
    echo "----------------------------------------"
done

echo "✅ 批次結構預測完成!"
```

## 🤖 機器學習工作流程

### 圖像分類模型訓練
```python
#!/usr/bin/env python3
# image_classification_pipeline.py - 圖像分類完整流程

import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import tensorflow as tf
from tensorflow.keras import layers, models
import seaborn as sns

def setup_environment():
    """設置環境和 GPU"""
    print("🔧 設置執行環境")
    
    # 檢查 GPU 可用性
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        print(f"✅ 發現 {len(gpus)} 個 GPU")
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    else:
        print("ℹ️  使用 CPU 模式")
    
    # 設置隨機種子
    tf.random.set_seed(42)
    np.random.seed(42)

def load_and_preprocess_data(data_dir, img_size=(224, 224), batch_size=32):
    """載入和預處理資料"""
    print("📊 載入和預處理資料")
    
    # 資料增強
    train_datagen = tf.keras.preprocessing.image.ImageDataGenerator(
        rescale=1./255,
        rotation_range=20,
        width_shift_range=0.2,
        height_shift_range=0.2,
        horizontal_flip=True,
        validation_split=0.2
    )
    
    # 訓練資料
    train_data = train_datagen.flow_from_directory(
        data_dir,
        target_size=img_size,
        batch_size=batch_size,
        class_mode='categorical',
        subset='training'
    )
    
    # 驗證資料
    val_data = train_datagen.flow_from_directory(
        data_dir,
        target_size=img_size,
        batch_size=batch_size,
        class_mode='categorical',
        subset='validation'
    )
    
    return train_data, val_data

def create_model(input_shape, num_classes):
    """建立 CNN 模型"""
    print("🏗️  建立 CNN 模型")
    
    model = models.Sequential([
        layers.Conv2D(32, (3, 3), activation='relu', input_shape=input_shape),
        layers.MaxPooling2D((2, 2)),
        layers.Conv2D(64, (3, 3), activation='relu'),
        layers.MaxPooling2D((2, 2)),
        layers.Conv2D(128, (3, 3), activation='relu'),
        layers.MaxPooling2D((2, 2)),
        layers.Conv2D(128, (3, 3), activation='relu'),
        layers.MaxPooling2D((2, 2)),
        layers.Flatten(),
        layers.Dropout(0.5),
        layers.Dense(512, activation='relu'),
        layers.Dense(num_classes, activation='softmax')
    ])
    
    model.compile(
        optimizer='adam',
        loss='categorical_crossentropy',
        metrics=['accuracy']
    )
    
    return model

def train_model(model, train_data, val_data, epochs=20):
    """訓練模型"""
    print("🚀 開始模型訓練")
    
    # 回調函數
    callbacks = [
        tf.keras.callbacks.EarlyStopping(
            monitor='val_loss',
            patience=5,
            restore_best_weights=True
        ),
        tf.keras.callbacks.ReduceLROnPlateau(
            monitor='val_loss',
            factor=0.2,
            patience=3,
            min_lr=0.0001
        ),
        tf.keras.callbacks.ModelCheckpoint(
            'results/models/best_model.h5',
            monitor='val_accuracy',
            save_best_only=True
        )
    ]
    
    # 訓練
    history = model.fit(
        train_data,
        validation_data=val_data,
        epochs=epochs,
        callbacks=callbacks,
        verbose=1
    )
    
    return history

def evaluate_and_visualize(model, val_data, class_names, output_dir):
    """評估模型並生成視覺化"""
    print("📊 評估模型效能")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 預測
    predictions = model.predict(val_data)
    y_pred = np.argmax(predictions, axis=1)
    y_true = val_data.classes
    
    # 分類報告
    report = classification_report(y_true, y_pred, target_names=class_names)
    print("\n分類報告:")
    print(report)
    
    # 儲存報告
    with open(f"{output_dir}/classification_report.txt", "w") as f:
        f.write(report)
    
    # 混淆矩陣
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=class_names, yticklabels=class_names)
    plt.title('混淆矩陣')
    plt.ylabel('真實標籤')
    plt.xlabel('預測標籤')
    plt.savefig(f"{output_dir}/confusion_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    return report

def plot_training_history(history, output_dir):
    """繪製訓練歷史"""
    print("📈 生成訓練歷史圖表")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # 準確率
    ax1.plot(history.history['accuracy'], label='訓練準確率')
    ax1.plot(history.history['val_accuracy'], label='驗證準確率')
    ax1.set_title('模型準確率')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('準確率')
    ax1.legend()
    
    # 損失
    ax2.plot(history.history['loss'], label='訓練損失')
    ax2.plot(history.history['val_loss'], label='驗證損失')
    ax2.set_title('模型損失')
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('損失')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/training_history.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """主要執行流程"""
    print("🤖 圖像分類機器學習流程")
    print("=" * 40)
    
    # 設定參數
    DATA_DIR = "data/images"  # 圖像資料目錄
    OUTPUT_DIR = "results/image_classification"
    IMG_SIZE = (224, 224)
    BATCH_SIZE = 32
    EPOCHS = 20
    
    # 檢查資料目錄
    if not os.path.exists(DATA_DIR):
        print(f"❌ 資料目錄不存在: {DATA_DIR}")
        print("請確保資料目錄結構如下:")
        print("data/images/")
        print("├── class1/")
        print("├── class2/")
        print("└── class3/")
        return
    
    # 建立輸出目錄
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(f"{OUTPUT_DIR}/models", exist_ok=True)
    
    # 執行流程
    setup_environment()
    
    train_data, val_data = load_and_preprocess_data(
        DATA_DIR, IMG_SIZE, BATCH_SIZE
    )
    
    # 取得類別資訊
    class_names = list(train_data.class_indices.keys())
    num_classes = len(class_names)
    print(f"類別數量: {num_classes}")
    print(f"類別名稱: {class_names}")
    
    # 建立和訓練模型
    model = create_model((*IMG_SIZE, 3), num_classes)
    print(model.summary())
    
    history = train_model(model, train_data, val_data, EPOCHS)
    
    # 評估和視覺化
    evaluate_and_visualize(model, val_data, class_names, OUTPUT_DIR)
    plot_training_history(history, OUTPUT_DIR)
    
    # 儲存最終模型
    model.save(f"{OUTPUT_DIR}/models/final_model.h5")
    
    print("✅ 圖像分類流程完成!")
    print(f"📂 結果位置: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
```

### 自然語言處理流程
```python
#!/usr/bin/env python3
# nlp_sentiment_analysis.py - 情感分析流程

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from transformers import TrainingArguments, Trainer
import torch
from torch.utils.data import Dataset
import os

class SentimentDataset(Dataset):
    """自訂資料集類別"""
    def __init__(self, texts, labels, tokenizer, max_length=512):
        self.texts = texts
        self.labels = labels
        self.tokenizer = tokenizer
        self.max_length = max_length
    
    def __len__(self):
        return len(self.texts)
    
    def __getitem__(self, idx):
        text = str(self.texts[idx])
        encoding = self.tokenizer(
            text,
            truncation=True,
            padding='max_length',
            max_length=self.max_length,
            return_tensors='pt'
        )
        
        return {
            'input_ids': encoding['input_ids'].flatten(),
            'attention_mask': encoding['attention_mask'].flatten(),
            'labels': torch.tensor(self.labels[idx], dtype=torch.long)
        }

def load_and_preprocess_data(data_file):
    """載入和預處理文本資料"""
    print("📊 載入和預處理文本資料")
    
    # 載入 CSV 資料 (假設有 'text' 和 'sentiment' 欄位)
    df = pd.read_csv(data_file)
    
    # 資料清理
    df = df.dropna(subset=['text', 'sentiment'])
    df['text'] = df['text'].astype(str)
    
    # 標籤編碼
    label_mapping = {'negative': 0, 'neutral': 1, 'positive': 2}
    df['label'] = df['sentiment'].map(label_mapping)
    
    # 移除未知標籤
    df = df.dropna(subset=['label'])
    df['label'] = df['label'].astype(int)
    
    print(f"資料量: {len(df)} 條")
    print("標籤分布:")
    print(df['sentiment'].value_counts())
    
    return df

def setup_model_and_tokenizer(model_name="bert-base-chinese"):
    """設置模型和分詞器"""
    print(f"🤖 載入預訓練模型: {model_name}")
    
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSequenceClassification.from_pretrained(
        model_name, 
        num_labels=3
    )
    
    return model, tokenizer

def train_sentiment_model(model, train_dataset, val_dataset, output_dir):
    """訓練情感分析模型"""
    print("🚀 開始模型微調")
    
    training_args = TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=3,
        per_device_train_batch_size=16,
        per_device_eval_batch_size=16,
        warmup_steps=500,
        weight_decay=0.01,
        logging_dir=f'{output_dir}/logs',
        logging_steps=100,
        evaluation_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="eval_loss",
    )
    
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
    )
    
    trainer.train()
    return trainer

def evaluate_model(trainer, test_dataset, output_dir):
    """評估模型"""
    print("📊 評估模型效能")
    
    # 預測
    predictions = trainer.predict(test_dataset)
    y_pred = np.argmax(predictions.predictions, axis=1)
    y_true = predictions.label_ids
    
    # 分類報告
    label_names = ['Negative', 'Neutral', 'Positive']
    report = classification_report(y_true, y_pred, target_names=label_names)
    print("\n分類報告:")
    print(report)
    
    # 儲存報告
    with open(f"{output_dir}/evaluation_report.txt", "w") as f:
        f.write(report)
    
    # 混淆矩陣
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=label_names, yticklabels=label_names)
    plt.title('情感分析混淆矩陣')
    plt.ylabel('真實標籤')
    plt.xlabel('預測標籤')
    plt.savefig(f"{output_dir}/confusion_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """主要執行流程"""
    print("🤖 自然語言處理 - 情感分析流程")
    print("=" * 40)
    
    # 設定參數
    DATA_FILE = "data/sentiment_data.csv"
    OUTPUT_DIR = "results/sentiment_analysis"
    MODEL_NAME = "bert-base-chinese"  # 或使用 "bert-base-uncased" 對於英文
    
    # 檢查資料檔案
    if not os.path.exists(DATA_FILE):
        print(f"❌ 資料檔案不存在: {DATA_FILE}")
        print("請準備包含 'text' 和 'sentiment' 欄位的 CSV 檔案")
        return
    
    # 建立輸出目錄
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 執行流程
    df = load_and_preprocess_data(DATA_FILE)
    
    # 分割資料
    train_df, test_df = train_test_split(df, test_size=0.2, random_state=42, stratify=df['label'])
    train_df, val_df = train_test_split(train_df, test_size=0.2, random_state=42, stratify=train_df['label'])
    
    print(f"訓練集: {len(train_df)} 條")
    print(f"驗證集: {len(val_df)} 條")
    print(f"測試集: {len(test_df)} 條")
    
    # 設置模型
    model, tokenizer = setup_model_and_tokenizer(MODEL_NAME)
    
    # 建立資料集
    train_dataset = SentimentDataset(
        train_df['text'].values, train_df['label'].values, tokenizer
    )
    val_dataset = SentimentDataset(
        val_df['text'].values, val_df['label'].values, tokenizer
    )
    test_dataset = SentimentDataset(
        test_df['text'].values, test_df['label'].values, tokenizer
    )
    
    # 訓練模型
    trainer = train_sentiment_model(model, train_dataset, val_dataset, OUTPUT_DIR)
    
    # 評估模型
    evaluate_model(trainer, test_dataset, OUTPUT_DIR)
    
    # 儲存模型
    trainer.save_model(f"{OUTPUT_DIR}/final_model")
    tokenizer.save_pretrained(f"{OUTPUT_DIR}/final_model")
    
    print("✅ 情感分析流程完成!")
    print(f"📂 結果位置: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
```

## 📊 資料科學分析流程

### 探索性資料分析 (EDA)
```python
#!/usr/bin/env python3
# exploratory_data_analysis.py - 探索性資料分析流程

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

def load_and_inspect_data(data_file):
    """載入和初步檢視資料"""
    print("📊 載入和初步檢視資料")
    
    # 載入資料
    df = pd.read_csv(data_file)
    
    print(f"資料形狀: {df.shape}")
    print(f"欄位數量: {df.shape[1]}")
    print(f"記錄數量: {df.shape[0]}")
    
    print("\n基本資訊:")
    print(df.info())
    
    print("\n前5筆記錄:")
    print(df.head())
    
    print("\n數值欄位統計:")
    print(df.describe())
    
    return df

def data_quality_analysis(df, output_dir):
    """資料品質分析"""
    print("🔍 資料品質分析")
    
    # 缺失值分析
    missing_data = df.isnull().sum()
    missing_percent = 100 * missing_data / len(df)
    
    missing_df = pd.DataFrame({
        '缺失數量': missing_data,
        '缺失百分比': missing_percent
    })
    missing_df = missing_df[missing_df['缺失數量'] > 0].sort_values('缺失數量', ascending=False)
    
    if len(missing_df) > 0:
        print("\n缺失值統計:")
        print(missing_df)
        
        # 缺失值視覺化
        plt.figure(figsize=(10, 6))
        missing_df['缺失百分比'].plot(kind='bar')
        plt.title('各欄位缺失值百分比')
        plt.ylabel('缺失百分比 (%)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/missing_values.png", dpi=300, bbox_inches='tight')
        plt.close()
    else:
        print("✅ 無缺失值")
    
    # 重複值檢查
    duplicates = df.duplicated().sum()
    print(f"\n重複記錄: {duplicates} 筆 ({100*duplicates/len(df):.2f}%)")
    
    # 資料類型分析
    print("\n資料類型分布:")
    print(df.dtypes.value_counts())
    
    return missing_df

def numerical_analysis(df, output_dir):
    """數值變數分析"""
    print("📈 數值變數分析")
    
    # 選擇數值欄位
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    if len(numeric_cols) == 0:
        print("無數值欄位")
        return
    
    print(f"數值欄位: {numeric_cols}")
    
    # 分布分析
    n_cols = min(3, len(numeric_cols))
    n_rows = (len(numeric_cols) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    for i, col in enumerate(numeric_cols):
        if i < len(axes):
            axes[i].hist(df[col].dropna(), bins=30, alpha=0.7)
            axes[i].set_title(f'{col} 分布')
            axes[i].set_xlabel(col)
            axes[i].set_ylabel('頻率')
    
    # 隱藏多餘的子圖
    for i in range(len(numeric_cols), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/numerical_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 相關性分析
    if len(numeric_cols) > 1:
        correlation_matrix = df[numeric_cols].corr()
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0,
                    square=True, fmt='.2f')
        plt.title('數值變數相關性矩陣')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/correlation_matrix.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 找出高相關性配對
        high_corr_pairs = []
        for i in range(len(correlation_matrix.columns)):
            for j in range(i+1, len(correlation_matrix.columns)):
                corr_val = correlation_matrix.iloc[i, j]
                if abs(corr_val) > 0.7:  # 高相關性閾值
                    high_corr_pairs.append((
                        correlation_matrix.columns[i],
                        correlation_matrix.columns[j],
                        corr_val
                    ))
        
        if high_corr_pairs:
            print("\n高相關性變數配對 (|r| > 0.7):")
            for var1, var2, corr in high_corr_pairs:
                print(f"  {var1} - {var2}: {corr:.3f}")

def categorical_analysis(df, output_dir):
    """類別變數分析"""
    print("📊 類別變數分析")
    
    # 選擇類別欄位
    categorical_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()
    
    if len(categorical_cols) == 0:
        print("無類別欄位")
        return
    
    print(f"類別欄位: {categorical_cols}")
    
    # 分析每個類別變數
    for col in categorical_cols[:4]:  # 限制最多4個欄位
        print(f"\n{col} 分布:")
        value_counts = df[col].value_counts()
        print(value_counts.head(10))
        
        # 長條圖
        plt.figure(figsize=(10, 6))
        if len(value_counts) <= 20:  # 類別不太多時顯示所有類別
            value_counts.plot(kind='bar')
            plt.xticks(rotation=45)
        else:  # 類別太多時只顯示前20個
            value_counts.head(20).plot(kind='bar')
            plt.xticks(rotation=45)
            plt.title(f'{col} 分布 (前20個類別)')
        
        plt.title(f'{col} 分布')
        plt.ylabel('頻率')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/categorical_{col.replace('/', '_')}.png", dpi=300, bbox_inches='tight')
        plt.close()

def outlier_analysis(df, output_dir):
    """異常值分析"""
    print("🎯 異常值分析")
    
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    
    if len(numeric_cols) == 0:
        print("無數值欄位進行異常值分析")
        return
    
    # 箱型圖
    n_cols = min(3, len(numeric_cols))
    n_rows = (len(numeric_cols) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    outlier_summary = {}
    
    for i, col in enumerate(numeric_cols):
        if i < len(axes):
            # 箱型圖
            axes[i].boxplot(df[col].dropna())
            axes[i].set_title(f'{col} 箱型圖')
            axes[i].set_ylabel(col)
            
            # IQR 方法檢測異常值
            Q1 = df[col].quantile(0.25)
            Q3 = df[col].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            
            outliers = df[(df[col] < lower_bound) | (df[col] > upper_bound)][col]
            outlier_summary[col] = {
                '異常值數量': len(outliers),
                '異常值比例': len(outliers) / len(df) * 100,
                '下界': lower_bound,
                '上界': upper_bound
            }
    
    # 隱藏多餘的子圖
    for i in range(len(numeric_cols), len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/outlier_boxplots.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 異常值摘要
    print("\n異常值摘要 (IQR 方法):")
    outlier_df = pd.DataFrame(outlier_summary).T
    print(outlier_df)
    
    return outlier_df

def generate_summary_report(df, missing_df, outlier_df, output_dir):
    """生成摘要報告"""
    print("📋 生成摘要報告")
    
    report = f"""
========================================
探索性資料分析報告
========================================
分析日期: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

資料概要:
- 總記錄數: {df.shape[0]:,}
- 總欄位數: {df.shape[1]}
- 數值欄位: {len(df.select_dtypes(include=[np.number]).columns)}
- 類別欄位: {len(df.select_dtypes(include=['object', 'category']).columns)}

資料品質:
- 完整記錄: {df.dropna().shape[0]:,} ({100*df.dropna().shape[0]/df.shape[0]:.1f}%)
- 有缺失值記錄: {df.isnull().any(axis=1).sum():,} ({100*df.isnull().any(axis=1).sum()/df.shape[0]:.1f}%)
- 重複記錄: {df.duplicated().sum():,} ({100*df.duplicated().sum()/df.shape[0]:.1f}%)

缺失值最嚴重的欄位:
"""
    
    if len(missing_df) > 0:
        for idx, row in missing_df.head(5).iterrows():
            report += f"- {idx}: {row['缺失數量']} ({row['缺失百分比']:.1f}%)\n"
    else:
        report += "- 無缺失值\n"
    
    report += "\n異常值最多的欄位:\n"
    if len(outlier_df) > 0:
        for idx, row in outlier_df.sort_values('異常值數量', ascending=False).head(5).iterrows():
            report += f"- {idx}: {row['異常值數量']} ({row['異常值比例']:.1f}%)\n"
    else:
        report += "- 無數值欄位進行異常值分析\n"
    
    report += f"""
建議後續分析:
1. 處理缺失值 (填補或移除)
2. 處理異常值 (轉換或移除)  
3. 特徵工程 (編碼、標準化)
4. 進一步的統計分析或建模

生成的視覺化檔案:
- missing_values.png: 缺失值分析
- numerical_distributions.png: 數值變數分布
- correlation_matrix.png: 相關性矩陣
- categorical_*.png: 類別變數分布
- outlier_boxplots.png: 異常值箱型圖
========================================
"""
    
    # 儲存報告
    with open(f"{output_dir}/eda_report.txt", "w", encoding='utf-8') as f:
        f.write(report)
    
    print("✅ 摘要報告已儲存")

def main():
    """主要執行流程"""
    print("📊 探索性資料分析流程")
    print("=" * 40)
    
    # 設定參數
    DATA_FILE = "data/dataset.csv"  # 修改為您的資料檔案路徑
    OUTPUT_DIR = "results/eda"
    
    # 檢查資料檔案
    if not os.path.exists(DATA_FILE):
        print(f"❌ 資料檔案不存在: {DATA_FILE}")
        print("請準備 CSV 格式的資料檔案")
        return
    
    # 建立輸出目錄
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 執行分析流程
    df = load_and_inspect_data(DATA_FILE)
    missing_df = data_quality_analysis(df, OUTPUT_DIR)
    numerical_analysis(df, OUTPUT_DIR)
    categorical_analysis(df, OUTPUT_DIR)
    outlier_df = outlier_analysis(df, OUTPUT_DIR)
    generate_summary_report(df, missing_df, outlier_df, OUTPUT_DIR)
    
    print("✅ 探索性資料分析完成!")
    print(f"📂 結果位置: {OUTPUT_DIR}")
    print(f"📋 查看報告: {OUTPUT_DIR}/eda_report.txt")

if __name__ == "__main__":
    main()
```

## 🚀 工作流程自動化

### 通用分析流程管理器
```bash
#!/bin/bash
# workflow_manager.sh - 通用分析流程管理器

set -e

WORKFLOW_DIR="workflows"
RESULTS_DIR="results"
LOG_DIR="logs"

# 顏色定義
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日誌函數
log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

# 建立目錄結構
setup_directories() {
    log_info "建立工作流程目錄結構"
    
    mkdir -p "$WORKFLOW_DIR"
    mkdir -p "$RESULTS_DIR"
    mkdir -p "$LOG_DIR"
    
    # 建立標準資料目錄
    mkdir -p data/{raw,processed,reference}
    mkdir -p scripts/{preprocessing,analysis,postprocessing}
    mkdir -p configs
    
    log_success "目錄結構建立完成"
}

# 列出可用工作流程
list_workflows() {
    echo "🔍 可用的工作流程:"
    echo "=================="
    
    echo "1. ngs          - NGS 基因組分析流程"
    echo "2. rnaseq       - RNA-seq 表達分析流程"  
    echo "3. structure    - 蛋白質結構預測流程"
    echo "4. ml-image     - 機器學習圖像分類流程"
    echo "5. ml-nlp       - 自然語言處理流程"
    echo "6. eda          - 探索性資料分析流程"
    echo "7. custom       - 自訂工作流程"
    echo ""
}

# 執行 NGS 工作流程
run_ngs_workflow() {
    local sample_name="$1"
    local fastq_dir="$2"
    local reference="$3"
    
    log_info "啟動 NGS 工作流程: $sample_name"
    
    # 檢查必要檔案
    if [ ! -d "$fastq_dir" ]; then
        log_error "FASTQ 目錄不存在: $fastq_dir"
        return 1
    fi
    
    if [ ! -f "$reference" ]; then
        log_error "參考基因組檔案不存在: $reference"
        return 1
    fi
    
    # 執行分析
    log_info "執行完整 NGS 分析流程..."
    bash workflows/wgs_analysis.sh "$sample_name" "$fastq_dir" "$reference" 2>&1 | tee "$LOG_DIR/ngs_${sample_name}_$(date +%Y%m%d_%H%M%S).log"
    
    log_success "NGS 工作流程完成: $sample_name"
}

# 執行結構預測工作流程
run_structure_workflow() {
    local protein_fasta="$1"
    local protein_name="$2"
    
    log_info "啟動蛋白質結構預測工作流程: $protein_name"
    
    if [ ! -f "$protein_fasta" ]; then
        log_error "蛋白質序列檔案不存在: $protein_fasta"
        return 1
    fi
    
    # 執行結構預測
    log_info "執行 ColabFold 結構預測..."
    bash workflows/protein_structure_prediction.sh "$protein_fasta" "$protein_name" 2>&1 | tee "$LOG_DIR/structure_${protein_name}_$(date +%Y%m%d_%H%M%S).log"
    
    log_success "結構預測工作流程完成: $protein_name"
}

# 執行機器學習工作流程
run_ml_workflow() {
    local workflow_type="$1"
    local data_dir="$2"
    
    log_info "啟動機器學習工作流程: $workflow_type"
    
    case "$workflow_type" in
        "image")
            log_info "執行圖像分類流程..."
            python3 workflows/image_classification_pipeline.py 2>&1 | tee "$LOG_DIR/ml_image_$(date +%Y%m%d_%H%M%S).log"
            ;;
        "nlp")
            log_info "執行自然語言處理流程..."
            python3 workflows/nlp_sentiment_analysis.py 2>&1 | tee "$LOG_DIR/ml_nlp_$(date +%Y%m%d_%H%M%S).log"
            ;;
        *)
            log_error "未知的機器學習工作流程類型: $workflow_type"
            return 1
            ;;
    esac
    
    log_success "機器學習工作流程完成: $workflow_type"
}

# 執行 EDA 工作流程
run_eda_workflow() {
    local data_file="$1"
    
    log_info "啟動探索性資料分析工作流程"
    
    if [ ! -f "$data_file" ]; then
        log_error "資料檔案不存在: $data_file"
        return 1
    fi
    
    # 執行 EDA
    log_info "執行探索性資料分析..."
    python3 workflows/exploratory_data_analysis.py 2>&1 | tee "$LOG_DIR/eda_$(date +%Y%m%d_%H%M%S).log"
    
    log_success "EDA 工作流程完成"
}

# 工作流程狀態檢查
check_workflow_status() {
    local workflow_id="$1"
    
    log_info "檢查工作流程狀態: $workflow_id"
    
    # 檢查結果目錄
    if [ -d "$RESULTS_DIR/$workflow_id" ]; then
        echo "📂 結果目錄: $RESULTS_DIR/$workflow_id"
        echo "📊 檔案數量: $(find "$RESULTS_DIR/$workflow_id" -type f | wc -l)"
        echo "💾 目錄大小: $(du -sh "$RESULTS_DIR/$workflow_id" | cut -f1)"
    else
        log_warning "結果目錄不存在: $RESULTS_DIR/$workflow_id"
    fi
    
    # 檢查日誌檔案
    local latest_log=$(ls -t "$LOG_DIR"/*"$workflow_id"* 2>/dev/null | head -1)
    if [ -n "$latest_log" ]; then
        echo "📋 最新日誌: $latest_log"
        echo "⏰ 日誌時間: $(stat -c %y "$latest_log")"
    else
        log_warning "找不到相關日誌檔案"
    fi
}

# 清理工作流程結果
cleanup_workflow() {
    local workflow_id="$1"
    local confirm="$2"
    
    if [ "$confirm" != "yes" ]; then
        echo "⚠️  這將刪除工作流程 '$workflow_id' 的所有結果和日誌"
        echo "如要確認，請執行: $0 cleanup $workflow_id yes"
        return 1
    fi
    
    log_warning "清理工作流程結果: $workflow_id"
    
    # 刪除結果目錄
    if [ -d "$RESULTS_DIR/$workflow_id" ]; then
        rm -rf "$RESULTS_DIR/$workflow_id"
        log_info "已刪除結果目錄: $RESULTS_DIR/$workflow_id"
    fi
    
    # 刪除相關日誌
    find "$LOG_DIR" -name "*$workflow_id*" -type f -delete
    log_info "已刪除相關日誌檔案"
    
    log_success "清理完成: $workflow_id"
}

# 主要命令處理
main() {
    case "$1" in
        "setup")
            setup_directories
            ;;
        "list")
            list_workflows
            ;;
        "ngs")
            if [ $# -ne 4 ]; then
                echo "用法: $0 ngs <sample_name> <fastq_dir> <reference_genome>"
                exit 1
            fi
            run_ngs_workflow "$2" "$3" "$4"
            ;;
        "structure")
            if [ $# -ne 3 ]; then
                echo "用法: $0 structure <protein.fasta> <protein_name>"
                exit 1
            fi
            run_structure_workflow "$2" "$3"
            ;;
        "ml")
            if [ $# -ne 3 ]; then
                echo "用法: $0 ml <image|nlp> <data_dir>"
                exit 1
            fi
            run_ml_workflow "$2" "$3"
            ;;
        "eda")
            if [ $# -ne 2 ]; then
                echo "用法: $0 eda <data_file.csv>"
                exit 1
            fi
            run_eda_workflow "$2"
            ;;
        "status")
            if [ $# -ne 2 ]; then
                echo "用法: $0 status <workflow_id>"
                exit 1
            fi
            check_workflow_status "$2"
            ;;
        "cleanup")
            if [ $# -lt 2 ]; then
                echo "用法: $0 cleanup <workflow_id> [yes]"
                exit 1
            fi
            cleanup_workflow "$2" "$3"
            ;;
        *)
            echo "🚀 Lab Container 工作流程管理器"
            echo "================================"
            echo ""
            echo "用法: $0 <命令> [參數]"
            echo ""
            echo "可用命令:"
            echo "  setup                           - 建立工作流程目錄結構"
            echo "  list                            - 列出可用工作流程"
            echo "  ngs <sample> <fastq_dir> <ref>  - 執行 NGS 分析流程"
            echo "  structure <fasta> <name>        - 執行蛋白質結構預測"
            echo "  ml <image|nlp> <data_dir>       - 執行機器學習流程"
            echo "  eda <data_file.csv>             - 執行探索性資料分析"
            echo "  status <workflow_id>            - 檢查工作流程狀態"
            echo "  cleanup <workflow_id> [yes]     - 清理工作流程結果"
            echo ""
            echo "範例:"
            echo "  $0 setup"
            echo "  $0 ngs sample01 data/fastq data/ref/hg38.fa"
            echo "  $0 structure protein.fasta my_protein"
            echo "  $0 ml image data/images"
            echo "  $0 eda data/dataset.csv"
            ;;
    esac
}

# 執行主函數
main "$@"
```

---

**這些工作流程範例涵蓋了 Lab Container 系統的主要應用領域**  
*從基因組學到機器學習，從蛋白質結構到資料科學，提供完整的分析流程模板！*