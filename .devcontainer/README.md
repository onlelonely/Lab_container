# Dev Container 彈性環境系統

## 📋 概述

這是一個強大的 Dev Container 彈性環境系統，支援從基本資料科學到專業生物資訊分析的多種工作流程。系統提供預設環境配置、動態套件管理、專用工具容器管理等功能。

## 🎯 主要功能

### 1. 🔧 環境 Profiles 系統
- **minimal** - 基本 Python/R 環境
- **datascience** - 資料科學標準環境 (pandas, tidyverse)
- **bioinformatics** - NGS 生物資訊環境 (Bioconductor, BioPython, FastQC, SAMtools)
- **ml** - 機器學習環境 (scikit-learn, tensorflow)
- **statistics** - 統計分析環境 (R focus)
- **full** - 完整功能環境
- **custom** - 自訂環境

### 2. 🐳 混合容器架構
```
📦 Main Dev Container (輕量級)         📦 專用工具容器群
├── 🐍 Python/R 分析環境              ├── 🔬 GATK (變異呼叫)
├── 📊 Jupyter Lab/RStudio             ├── 🧬 ColabFold (結構預測)
├── 🧬 輕量級生物資訊工具              └── 💊 AutoDock (分子對接)
│   ├── FastQC, SAMtools               
│   ├── BCFtools, BEDtools             
│   └── BioPython, pysam               
└── 📝 資料分析腳本                   
```

### 3. 📦 動態套件管理
- 運行時新增/移除 Python、R、Conda 套件
- 自動記錄安裝歷史
- 環境備份和還原功能

### 4. 🔧 配置範本系統
- 自動生成 devcontainer.json, docker-compose.yml
- 根據 profile 自動調整 VS Code 擴展和設定
- 支援條件配置和變數替換

## 📁 目錄結構

```
.devcontainer/
├── devcontainer.json              # 自動生成的容器配置
├── docker-compose.yml             # 自動生成的服務編排
├── docker-compose.tools.yml       # 專用工具容器定義
├── Dockerfile.generated           # 自動生成的 Dockerfile
├── README.md                      # 本文檔
├── README-profiles.md             # 詳細使用指南
├── profiles/                      # 環境設定檔
│   ├── minimal.env
│   ├── datascience.env
│   ├── bioinformatics.env
│   ├── ml.env
│   ├── statistics.env
│   ├── full.env
│   └── custom.env                 # 自動生成
├── templates/                     # 配置範本
│   ├── devcontainer.json.template
│   └── docker-compose.yml.template
├── scripts/
│   ├── manage/
│   │   ├── setup-environment.sh   # 環境選擇器
│   │   ├── devcontainer-status.sh # 狀態監控
│   │   ├── quick-start.sh         # 快速啟動
│   │   └── tool-manager.sh        # 專用工具管理
│   └── utils/
│       ├── package-manager.sh     # 動態套件管理
│       ├── template-generator.sh  # 配置生成器
│       └── post-create.sh         # 容器建立後腳本
├── configs/                       # 自動生成的配置
│   ├── python/
│   │   └── requirements-*.txt
│   ├── r/
│   │   └── *-packages.R
│   └── conda/
│       └── environment-*.yml
└── tests/                         # 測試套件
    ├── test-installation.sh
    └── test-environment.sh
```

## 🚀 快速開始

### 初次設定
```bash
# 執行環境選擇器
bash .devcontainer/scripts/manage/setup-environment.sh

# 或使用快速啟動
bash .devcontainer/scripts/manage/quick-start.sh
```

### 選擇環境
```bash
# 切換到生物資訊環境
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics

# 切換到資料科學環境
bash .devcontainer/scripts/manage/setup-environment.sh datascience
```

## 🔧 主要工具

### 環境管理
```bash
# 環境設定
setup-env                       # 重新設定環境
status                          # 檢查環境狀態
quick-start                     # 快速啟動

# 配置生成
bash .devcontainer/scripts/utils/template-generator.sh
```

### 套件管理
```bash
# 動態新增套件
add-package python requests     # 新增 Python 套件
add-package r forecast          # 新增 R 套件
add-package conda nodejs        # 新增 Conda 套件

# 套件管理
list-packages                   # 列出套件
backup-env                      # 備份環境
```

### 專用工具管理
```bash
# 重量級生物資訊工具
bash .devcontainer/scripts/manage/tool-manager.sh gatk --version
bash .devcontainer/scripts/manage/tool-manager.sh colabfold protein.fasta results/
bash .devcontainer/scripts/manage/tool-manager.sh autodock --receptor protein.pdbqt

# 工具容器管理
bash .devcontainer/scripts/manage/tool-manager.sh manage start
bash .devcontainer/scripts/manage/tool-manager.sh manage status
bash .devcontainer/scripts/manage/tool-manager.sh check all
```

## 🧬 生物資訊工作流程

### NGS 變異分析
```bash
# 1. 品質控制
fastqc data/fastq/*.fastq.gz -o results/qc/

# 2. 序列比對
bwa mem reference.fa sample.fastq.gz | samtools sort -o sample.bam

# 3. 變異呼叫 (專用容器)
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I sample.bam -O variants.vcf -R reference.fa

# 4. 變異過濾
bcftools filter variants.vcf > filtered_variants.vcf
```

### 蛋白質結構預測
```bash
# 使用 ColabFold (支援 GPU 加速)
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta results/structure/ --num-models 5
```

### 分子對接
```bash
# 使用 AutoDock Vina
bash .devcontainer/scripts/manage/tool-manager.sh autodock \
    --receptor protein.pdbqt --ligand ligand.pdbqt --out docking_result.pdbqt
```

## 📊 狀態監控

### 環境狀態檢查
```bash
# 完整狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh

# 快速狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh quick
```

## 🎯 特色功能

### 安全性
- 非根用戶 (vscode) 配置
- 適當的檔案權限和擁有權
- 安全的容器配置

### 效能優化
- 模組化配置檔案
- 重試邏輯和指數退避
- 完整的錯誤處理和日誌記錄
- 條件式工具載入

### 開發工具整合
- Python 3.11 與資料科學套件庫
- R 4.3 與 Bioconductor 支援
- Jupyter Lab 與多核心支援
- VS Code 擴展自動配置

## 📚 詳細文檔

更多詳細資訊請參考：
- **README-profiles.md** - 完整使用指南和範例
- **improve.md** - 改進建議和技術細節
- 各腳本內建的 help 功能

## 🔗 相關連結

- [VS Code Dev Containers](https://code.visualstudio.com/docs/remote/containers)
- [Docker Compose](https://docs.docker.com/compose/)
- [Bioconductor](https://bioconductor.org/)
- [GATK](https://gatk.broadinstitute.org/)
- [ColabFold](https://colab.research.google.com/github/deepmind/alphafold)
- [AutoDock Vina](https://autodock-vina.readthedocs.io/)
