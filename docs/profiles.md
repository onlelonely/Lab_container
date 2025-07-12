# Profile 系統說明

**文檔版本**: v2.0  
**系統組件**: Environment Profile System  
**適用讀者**: 所有使用者

## 🎯 Profile 系統概述

Profile 系統是 Lab Container 的核心功能，提供**預設環境模板**，讓使用者能夠快速設置特定用途的分析環境。

### 設計理念
- **專業化**: 每個 Profile 針對特定分析領域優化
- **標準化**: 提供業界標準的套件組合
- **可客製化**: 支援在 Profile 基礎上進一步客製
- **一致性**: 確保團隊間環境一致性

## 📋 可用 Profile 清單

### 1. 🟢 Minimal Profile
**適用對象**: 初學者、資源受限環境

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh minimal
```

**特色**:
- 最小資源需求 (< 2GB RAM)
- 基礎 Python 3.11 + R 4.3
- 快速啟動 (< 5 分鐘)
- 適合學習和簡單分析

**包含套件**:
```python
# Python 套件
numpy>=1.24.0
pandas>=2.0.0  
matplotlib>=3.7.0

# R 套件
base, utils, stats, graphics
```

**使用案例**:
- Python/R 語言學習
- 簡單資料處理
- 環境測試
- CI/CD 輕量化建構

### 2. 📊 Data Science Profile
**適用對象**: 資料科學家、分析師

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh datascience
```

**特色**:
- 標準資料科學工具鏈
- Jupyter Lab 完整整合
- 視覺化工具齊全
- 平衡效能與功能

**包含套件**:
```python
# Python 套件
pandas>=2.0.0        # 資料處理
numpy>=1.24.0         # 數值計算
matplotlib>=3.7.0     # 基礎繪圖
seaborn>=0.12.0       # 統計視覺化
plotly>=5.17.0        # 互動式圖表
scikit-learn>=1.3.0   # 機器學習
jupyterlab>=4.0.0     # 互動式開發

# R 套件  
tidyverse             # R 資料科學工具鏈
ggplot2               # 進階視覺化
dplyr                 # 資料操作
readr                 # 資料讀取
```

**使用案例**:
- 探索性資料分析 (EDA)
- 商業智慧報告
- 資料視覺化
- 統計分析

### 3. 🧬 Bioinformatics Profile  
**適用對象**: 生物資訊學研究者、NGS 分析

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics
```

**特色**:
- NGS 資料分析工具鏈
- 生物資料庫整合
- 基因組學工具
- Bioconductor 生態系統

**包含套件**:
```python
# Python 生物資訊套件
biopython>=1.81       # 生物序列分析
pysam>=0.22.0         # SAM/BAM 檔案處理
pyvcf>=0.6.8          # VCF 檔案處理

# 輕量級命令列工具 (透過 Conda)
fastqc                # 序列品質控制
samtools              # SAM/BAM 操作工具
bcftools              # VCF 檔案工具
bedtools              # 基因組區間工具

# R Bioconductor 套件
Biostrings            # 生物序列操作
GenomicRanges         # 基因組區間分析
VariantAnnotation     # 變異註釋
```

**專用工具容器** (可選):
```bash
# GATK 變異呼叫工具
bash .devcontainer/scripts/manage/tool-manager.sh gatk --version

# 更多重量級工具請參考工具管理文檔
```

**使用案例**:
- NGS 資料品質控制
- 基因組比對分析
- 變異檢測與註釋
- 基因表達分析

### 4. 🤖 Machine Learning Profile
**適用對象**: 機器學習工程師、AI 研究者

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh ml
```

**特色**:
- 深度學習框架齊全
- GPU 支援 (如可用)
- 模型開發工具鏈
- 實驗追蹤整合

**包含套件**:
```python
# 機器學習框架
scikit-learn>=1.3.0   # 傳統機器學習
tensorflow>=2.13.0    # 深度學習框架
torch>=2.0.0          # PyTorch 深度學習
transformers>=4.30.0  # 預訓練模型

# 資料處理
pandas>=2.0.0         # 資料處理
numpy>=1.24.0         # 數值計算
scipy>=1.10.0         # 科學計算

# 視覺化與監控
matplotlib>=3.7.0     # 繪圖
seaborn>=0.12.0       # 統計視覺化  
tensorboard>=2.13.0   # 實驗監控

# R 機器學習套件
caret                 # 分類回歸訓練
randomForest          # 隨機森林
nnet                  # 神經網路
```

**GPU 支援**:
```bash
# 檢查 GPU 可用性
nvidia-smi

# TensorFlow GPU 測試
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# PyTorch GPU 測試  
python -c "import torch; print(torch.cuda.is_available())"
```

**使用案例**:
- 深度學習模型訓練
- 電腦視覺專案
- 自然語言處理
- 推薦系統開發

### 5. 📈 Statistics Profile
**適用對象**: 統計學家、學術研究者

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh statistics
```

**特色**:
- R 語言為主的統計環境
- 進階統計套件
- 報告生成工具
- 學術出版支援

**包含套件**:
```r
# R 統計套件
stats                 # 基礎統計
MASS                  # 現代應用統計
survival              # 生存分析
forecast              # 時間序列預測
lme4                  # 線性混合效應模型
ggplot2               # 統計圖形
knitr                 # 動態報告
rmarkdown             # R Markdown 文檔

# Python 補充套件
scipy.stats           # 統計函數
statsmodels           # 統計模型
pingouin              # 統計檢定
```

**專用功能**:
- R Markdown 報告生成
- 學術圖表製作
- 統計檢定與建模
- 實驗設計分析

**使用案例**:
- 學術研究統計分析
- 臨床試驗分析
- 品質控制統計
- 統計教學

### 6. 🚀 Full Profile
**適用對象**: 進階使用者、多領域分析

```bash
# 啟動方式
bash .devcontainer/scripts/manage/setup-environment.sh full
```

**特色**:
- 包含所有上述 Profile 的功能
- 最大功能完整性
- 適合複雜多步驟分析
- 較大資源需求 (> 8GB RAM)

**包含所有套件**:
- Data Science Profile 所有套件
- Bioinformatics Profile 所有工具
- Machine Learning Profile 所有框架
- Statistics Profile 所有統計套件

**額外功能**:
- RStudio Server 整合
- 所有專用工具容器可用
- 完整的可視化能力
- 多語言混合分析支援

**使用案例**:
- 大型跨領域專案
- 教學與示範環境
- 全功能開發環境
- 生產環境部署

## 🛠 Profile 管理操作

### 查看當前 Profile
```bash
# 檢查當前啟用的 Profile
cat .devcontainer/.current-profile

# 詳細狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh
```

### 切換 Profile
```bash
# 重新選擇 Profile
bash .devcontainer/scripts/manage/setup-environment.sh

# 或直接指定 Profile
bash .devcontainer/scripts/manage/setup-environment.sh datascience
```

### Profile 客製化
```bash
# 在現有 Profile 基礎上安裝額外套件
bash .devcontainer/scripts/utils/package-manager.sh add python seaborn
bash .devcontainer/scripts/utils/package-manager.sh add r forecast

# 備份客製化環境
bash .devcontainer/scripts/utils/package-manager.sh backup
```

## 📊 Profile 比較表

| 功能特色 | Minimal | Data Science | Bioinformatics | ML | Statistics | Full |
|---------|---------|--------------|----------------|----|------------|------|
| **啟動時間** | < 5 min | < 10 min | < 15 min | < 20 min | < 10 min | < 25 min |
| **記憶體需求** | < 2GB | < 4GB | < 6GB | < 8GB | < 4GB | > 8GB |
| **磁碟空間** | < 3GB | < 5GB | < 8GB | < 12GB | < 5GB | > 15GB |
| **Python** | ✅ 基礎 | ✅ 完整 | ✅ 生物資訊 | ✅ ML/DL | ✅ 統計 | ✅ 全部 |
| **R** | ✅ 基礎 | ✅ tidyverse | ✅ Bioconductor | ✅ ML 套件 | ✅ 完整統計 | ✅ 全部 |
| **Jupyter** | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **RStudio** | ❌ | ❌ | ❌ | ❌ | ✅ | ✅ |
| **生物資訊工具** | ❌ | ❌ | ✅ | ❌ | ❌ | ✅ |
| **深度學習** | ❌ | 基礎 | ❌ | ✅ | ❌ | ✅ |
| **GPU 支援** | ❌ | ❌ | ❌ | ✅ | ❌ | ✅ |

## 🎯 選擇指南

### 新手使用者
```bash
# 推薦順序
1. minimal      # 學習基礎
2. datascience  # 開始分析
3. 其他 profile # 根據需求
```

### 專業使用者
```bash
# 根據領域直接選擇
資料科學 → datascience
生物資訊 → bioinformatics  
機器學習 → ml
統計分析 → statistics
多領域 → full
```

### 資源受限環境
```bash
# 優先考慮
minimal > datascience > statistics > bioinformatics > ml > full
```

## 🔧 自訂 Profile 開發

### 建立新 Profile
```bash
# 1. 複製現有 Profile
cp .devcontainer/profiles/datascience.env .devcontainer/profiles/myprofile.env

# 2. 編輯配置
vim .devcontainer/profiles/myprofile.env

# 3. 建立套件配置檔案
touch .devcontainer/configs/python/requirements-myprofile.txt
touch .devcontainer/configs/r/myprofile-packages.R
touch .devcontainer/configs/conda/environment-myprofile.yml

# 4. 測試新 Profile
bash .devcontainer/scripts/manage/setup-environment.sh myprofile
```

### Profile 配置範例
```bash
# myprofile.env
PROFILE_NAME="My Custom Profile"
PROFILE_DESCRIPTION="自訂分析環境"

# 套件定義
PYTHON_PACKAGES="pandas numpy custom-package"
R_PACKAGES="tidyverse custom-r-package"

# 功能開關
ENABLE_JUPYTER=true
ENABLE_RSTUDIO=false
ENABLE_BIOINFORMATICS=false

# 環境變數
CUSTOM_VAR="custom_value"
```

## ⚠️ 常見問題

### Profile 切換後套件消失
```bash
# 原因: 不同 Profile 有不同的套件配置
# 解決: 在新 Profile 中重新安裝需要的套件
bash .devcontainer/scripts/utils/package-manager.sh add python <package_name>
```

### 記憶體不足錯誤
```bash
# 原因: 選擇的 Profile 超出系統資源
# 解決: 切換到更輕量的 Profile
bash .devcontainer/scripts/manage/setup-environment.sh minimal
```

### Profile 載入失敗
```bash
# 檢查 Profile 檔案是否存在
ls .devcontainer/profiles/

# 檢查語法錯誤
bash -n .devcontainer/profiles/your-profile.env
```

---

**Profile 系統讓環境管理變得簡單而強大**  
*選擇適合的 Profile，立即開始您的分析工作！*