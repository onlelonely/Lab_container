# Dev Container 彈性環境系統

## 📋 概述

Dev Container 彈性環境系統讓您可以輕鬆選擇和切換不同的分析環境配置，支援從基本的資料科學到專業的生物資訊分析等多種工作流程。

## 🎯 主要功能

### 1. 預設環境 Profiles
- **minimal** - 基本 Python/R 環境
- **datascience** - 資料科學標準環境 (pandas, tidyverse)
- **bioinformatics** - NGS 生物資訊環境 (Bioconductor, BioPython, FastQC, SAMtools)
- **ml** - 機器學習環境 (scikit-learn, tensorflow)
- **statistics** - 統計分析環境 (R focus)
- **full** - 完整功能環境
- **custom** - 自訂環境

### 2. 動態套件管理
- 運行時新增/移除 Python、R、Conda 套件
- 自動記錄安裝歷史
- 環境備份和還原功能

### 3. 智能狀態監控
- 即時環境狀態檢查
- 套件清單和版本資訊
- 系統資源監控

## 🚀 快速開始

### 初次設定
```bash
# 執行環境選擇器
bash .devcontainer/scripts/manage/setup-environment.sh

# 或使用快速啟動
bash .devcontainer/scripts/manage/quick-start.sh
```

### 切換環境
```bash
# 切換到生物資訊環境
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics

# 切換到資料科學環境
bash .devcontainer/scripts/manage/setup-environment.sh datascience
```

## 📦 套件管理

### 新增套件
```bash
# 新增 Python 套件
bash .devcontainer/scripts/utils/package-manager.sh add python requests

# 新增 R 套件
bash .devcontainer/scripts/utils/package-manager.sh add r forecast

# 新增 Conda 套件
bash .devcontainer/scripts/utils/package-manager.sh add conda nodejs
```

### 查看套件
```bash
# 列出所有已安裝套件
bash .devcontainer/scripts/utils/package-manager.sh list

# 搜尋套件
bash .devcontainer/scripts/utils/package-manager.sh search tensorflow
```

### 環境備份
```bash
# 備份當前環境
bash .devcontainer/scripts/utils/package-manager.sh backup

# 還原環境
bash .devcontainer/scripts/utils/package-manager.sh restore backups/20250712_015234
```

## 📊 狀態檢查

### 完整狀態檢查
```bash
bash .devcontainer/scripts/manage/devcontainer-status.sh
```

### 快速狀態檢查
```bash
bash .devcontainer/scripts/manage/devcontainer-status.sh quick
```

## 🧬 生物資訊工具

當選擇 `bioinformatics` 或 `full` 環境時，系統會自動整合以下工具：

### 輕量級工具 (Container 內)
- **FastQC** - 序列品質控制
- **SAMtools** - SAM/BAM 檔案操作
- **BCFtools** - VCF 檔案處理
- **BEDtools** - 基因組區間操作

### Python 生物資訊套件
- **BioPython** - 生物序列分析
- **pysam** - SAM/BAM 檔案處理
- **pyvcf** - VCF 檔案解析

### R 生物資訊套件
- **Biostrings** - 生物序列操作
- **GenomicRanges** - 基因組區間分析
- **VariantAnnotation** - 變異註釋

## 📁 目錄結構

```
.devcontainer/
├── profiles/                    # 環境設定檔
│   ├── minimal.env
│   ├── datascience.env
│   ├── bioinformatics.env
│   ├── ml.env
│   ├── statistics.env
│   ├── full.env
│   └── custom.env              # 自動生成
├── scripts/
│   ├── manage/
│   │   ├── setup-environment.sh    # 環境選擇器
│   │   ├── devcontainer-status.sh  # 狀態檢查
│   │   └── quick-start.sh          # 快速啟動
│   └── utils/
│       └── package-manager.sh      # 套件管理工具
├── configs/                     # 自動生成的配置
│   ├── python/
│   │   └── requirements-*.txt
│   ├── r/
│   │   └── *-packages.R
│   └── conda/
│       └── environment-*.yml
└── .current-profile            # 當前環境記錄
```

## 🛠️ 自訂環境

您可以建立自己的環境設定檔：

1. 複製現有的 profile 檔案
2. 修改套件清單和設定
3. 使用環境選擇器載入

範例：
```bash
# 複製並修改
cp .devcontainer/profiles/datascience.env .devcontainer/profiles/myenv.env
# 編輯 myenv.env
# 載入自訂環境
bash .devcontainer/scripts/manage/setup-environment.sh myenv
```

## 🔧 快捷命令

設定完成後，可使用以下快捷命令（需重新啟動 shell）：

```bash
# 環境管理
setup-env          # 重新設定環境
status             # 檢查環境狀態
quick-start        # 快速啟動
qs                 # 快速狀態檢查

# 套件管理
add-package        # 新增套件
list-packages      # 列出套件
backup-env         # 備份環境
```

## 💡 使用技巧

### 1. 工作流程建議
1. 根據專案需求選擇適當的 profile
2. 使用狀態檢查確認環境正常
3. 需要額外套件時使用動態新增
4. 定期備份環境設定

### 2. 效能優化
- 選擇最小化的環境配置
- 定期清理不需要的套件
- 使用 conda 安裝大型套件

### 3. 故障排除
```bash
# 檢查環境狀態
bash .devcontainer/scripts/manage/devcontainer-status.sh

# 重新生成配置
bash .devcontainer/scripts/manage/setup-environment.sh <profile-name>

# 查看日誌
cat /tmp/jupyter.log  # Jupyter 日誌
```

## 📚 範例使用場景

### 資料科學專案
```bash
# 設定環境
bash .devcontainer/scripts/manage/setup-environment.sh datascience

# 新增特定套件
bash .devcontainer/scripts/utils/package-manager.sh add python plotly
bash .devcontainer/scripts/utils/package-manager.sh add r shiny
```

### NGS 分析專案
```bash
# 設定生物資訊環境
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics

# 檢查工具可用性
fastqc --version
samtools --version

# 執行品質控制
fastqc sample.fastq.gz -o results/
```

### 機器學習專案
```bash
# 設定 ML 環境
bash .devcontainer/scripts/manage/setup-environment.sh ml

# 新增深度學習套件
bash .devcontainer/scripts/utils/package-manager.sh add python pytorch-lightning
```

## 🎯 未來發展

規劃中的功能：
- 專用容器支援 (GATK, ColabFold, AutoDock Vina)
- GPU 加速環境
- 工作流程自動化
- 雲端環境同步

---

**需要協助？** 執行 `bash .devcontainer/scripts/manage/setup-environment.sh` 檢視完整選項，或查看 `improve.md` 了解詳細技術資訊。