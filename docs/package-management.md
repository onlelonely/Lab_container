# 動態套件管理指南

**文檔版本**: v2.0  
**核心組件**: Package Manager System  
**適用讀者**: 所有使用者

## 🚀 套件管理系統概述

Lab Container 提供強大的**動態套件管理系統**，支援在運行時安裝、移除、更新套件，無需重建容器。

### 核心特色
- **多語言支援**: Python、R、Conda 套件統一管理
- **運行時安裝**: 無需重啟容器即可安裝新套件
- **智慧相依性**: 自動處理套件相依性
- **環境備份**: 完整的環境狀態備份與還原
- **衝突檢測**: 自動偵測並解決套件衝突

## 📦 套件管理器使用方法

### 基本命令格式
```bash
# 套件管理器位置
bash .devcontainer/scripts/utils/package-manager.sh <操作> <套件類型> <套件名稱>
```

### 快速開始
```bash
# 檢視說明
bash .devcontainer/scripts/utils/package-manager.sh help

# 列出已安裝套件
bash .devcontainer/scripts/utils/package-manager.sh list

# 搜尋套件
bash .devcontainer/scripts/utils/package-manager.sh search tensorflow
```

## 🐍 Python 套件管理

### 安裝 Python 套件
```bash
# 基本安裝
bash .devcontainer/scripts/utils/package-manager.sh add python requests

# 指定版本安裝
bash .devcontainer/scripts/utils/package-manager.sh add python "pandas>=2.0.0,<2.2.0"

# 從 GitHub 安裝
bash .devcontainer/scripts/utils/package-manager.sh add python git+https://github.com/user/repo.git

# 批次安裝
bash .devcontainer/scripts/utils/package-manager.sh add python "requests beautifulsoup4 lxml"
```

### 移除 Python 套件
```bash
# 移除單一套件
bash .devcontainer/scripts/utils/package-manager.sh remove python requests

# 移除多個套件
bash .devcontainer/scripts/utils/package-manager.sh remove python "requests beautifulsoup4"
```

### Python 套件進階操作
```bash
# 更新套件到最新版本
python -m pip install --upgrade pandas

# 檢查套件資訊
python -m pip show pandas

# 列出過時套件
python -m pip list --outdated

# 安裝開發版本
bash .devcontainer/scripts/utils/package-manager.sh add python "package_name --pre"
```

### 常用 Python 套件範例
```bash
# 資料科學套件
bash .devcontainer/scripts/utils/package-manager.sh add python "pandas numpy matplotlib seaborn"

# 機器學習套件
bash .devcontainer/scripts/utils/package-manager.sh add python "scikit-learn tensorflow torch"

# 網頁爬蟲套件
bash .devcontainer/scripts/utils/package-manager.sh add python "requests beautifulsoup4 selenium"

# 生物資訊套件
bash .devcontainer/scripts/utils/package-manager.sh add python "biopython pysam pyvcf"
```

## 📊 R 套件管理

### 安裝 R 套件
```bash
# 從 CRAN 安裝
bash .devcontainer/scripts/utils/package-manager.sh add r ggplot2

# 從 Bioconductor 安裝
bash .devcontainer/scripts/utils/package-manager.sh add r Biostrings

# 從 GitHub 安裝
bash .devcontainer/scripts/utils/package-manager.sh add r "devtools::install_github('user/repo')"

# 批次安裝
bash .devcontainer/scripts/utils/package-manager.sh add r "ggplot2 dplyr tidyr"
```

### R 套件進階操作
```bash
# 在 R 控制台中直接管理
R
> install.packages("forecast")
> BiocManager::install("GenomicRanges")
> devtools::install_github("rstudio/shiny")
> quit()

# 檢查 R 套件資訊
R -e "packageVersion('ggplot2')"
R -e "installed.packages()[,c('Package','Version')]"
```

### 常用 R 套件範例
```bash
# Tidyverse 生態系統
bash .devcontainer/scripts/utils/package-manager.sh add r "tidyverse dplyr ggplot2 readr"

# 統計分析套件
bash .devcontainer/scripts/utils/package-manager.sh add r "forecast lme4 survival"

# 生物資訊套件 (需要 Bioconductor)
R -e "BiocManager::install(c('Biostrings', 'GenomicRanges', 'VariantAnnotation'))"

# 機器學習套件
bash .devcontainer/scripts/utils/package-manager.sh add r "caret randomForest nnet"
```

## 🔬 Conda 套件管理

### 安裝 Conda 套件
```bash
# 基本安裝
bash .devcontainer/scripts/utils/package-manager.sh add conda numpy

# 從特定 channel 安裝
bash .devcontainer/scripts/utils/package-manager.sh add conda "bioconda::samtools"

# 指定版本安裝
bash .devcontainer/scripts/utils/package-manager.sh add conda "python=3.11"
```

### Conda 環境管理
```bash
# 列出 conda 套件
conda list

# 更新所有套件
conda update --all

# 檢查 channel 優先級
conda config --show channels

# 搜尋套件
conda search tensorflow
```

### 常用 Conda 套件範例
```bash
# 生物資訊工具
bash .devcontainer/scripts/utils/package-manager.sh add conda "bioconda::fastqc bioconda::samtools"

# 科學計算套件
bash .devcontainer/scripts/utils/package-manager.sh add conda "numpy scipy matplotlib"

# 系統工具
bash .devcontainer/scripts/utils/package-manager.sh add conda "git wget curl"
```

## 🔍 套件搜尋與資訊查詢

### 搜尋功能
```bash
# 搜尋 Python 套件
bash .devcontainer/scripts/utils/package-manager.sh search python tensorflow

# 搜尋 R 套件
bash .devcontainer/scripts/utils/package-manager.sh search r ggplot

# 搜尋 Conda 套件
bash .devcontainer/scripts/utils/package-manager.sh search conda samtools
```

### 套件資訊查詢
```bash
# Python 套件詳細資訊
python -m pip show pandas

# R 套件資訊
R -e "help(package='ggplot2')"

# Conda 套件資訊
conda info pandas
```

### 相依性分析
```bash
# Python 套件相依性
python -m pip show pandas | grep "Requires"

# R 套件相依性
R -e "tools::package_dependencies('ggplot2')"

# Conda 套件相依性
conda info --json pandas | grep depends
```

## 📋 套件清單管理

### 列出已安裝套件
```bash
# 列出所有套件
bash .devcontainer/scripts/utils/package-manager.sh list

# 分別列出各類型套件
bash .devcontainer/scripts/utils/package-manager.sh list python
bash .devcontainer/scripts/utils/package-manager.sh list r
bash .devcontainer/scripts/utils/package-manager.sh list conda
```

### 匯出套件清單
```bash
# 匯出 Python 套件清單
python -m pip freeze > requirements.txt

# 匯出 R 套件清單
R -e "write.csv(installed.packages()[,c('Package','Version')], 'r-packages.csv')"

# 匯出 Conda 環境
conda env export > environment.yml
```

### 從清單安裝套件
```bash
# 從 requirements.txt 安裝
python -m pip install -r requirements.txt

# 從 R 套件清單安裝
R -e "packages <- read.csv('r-packages.csv'); install.packages(packages$Package)"

# 從 Conda 環境檔案安裝
conda env update --file environment.yml
```

## 💾 環境備份與還原

### 自動備份
```bash
# 建立完整環境備份
bash .devcontainer/scripts/utils/package-manager.sh backup

# 備份將儲存在 backups/YYYYMMDD_HHMMSS/ 目錄
# 包含：
# - Python 套件清單 (requirements.txt)
# - R 套件清單 (r-packages.csv)  
# - Conda 環境檔案 (environment.yml)
# - Profile 設定檔案
```

### 還原環境
```bash
# 列出可用備份
ls backups/

# 還原特定備份
bash .devcontainer/scripts/utils/package-manager.sh restore backups/20250712_143000

# 還原最新備份
bash .devcontainer/scripts/utils/package-manager.sh restore latest
```

### 定期備份建議
```bash
# 建議在以下時機進行備份：
1. 重大套件安裝前
2. 專案重要節點
3. 系統升級前
4. 定期例行備份 (每週/每月)

# 可設定自動備份 (crontab)
# 每週日備份一次
0 2 * * 0 cd /path/to/project && bash .devcontainer/scripts/utils/package-manager.sh backup
```

## ⚡ 效能優化

### 快取管理
```bash
# 清理 pip 快取
python -m pip cache purge

# 清理 conda 快取
conda clean --all

# 檢查快取大小
du -sh ~/.cache/pip ~/.conda/pkgs
```

### 平行安裝
```bash
# Python 平行安裝 (使用 --parallel)
python -m pip install --upgrade pip setuptools wheel
python -m pip install package1 package2 package3

# Conda 平行安裝
conda install -y package1 package2 package3
```

### 離線安裝
```bash
# 下載套件供離線安裝
python -m pip download pandas -d offline_packages/

# 離線安裝
python -m pip install --no-index --find-links offline_packages/ pandas
```

## 🚨 常見問題與解決方案

### 套件衝突
```bash
# 問題: 套件版本衝突
# 解決方案 1: 使用虛擬環境
python -m venv myenv
source myenv/bin/activate
pip install conflicting_package

# 解決方案 2: 強制重新安裝
pip install --force-reinstall package_name

# 解決方案 3: 檢查相依性
pip check
```

### 安裝失敗
```bash
# 問題: 編譯錯誤
# 解決方案: 使用 conda 或預編譯版本
conda install problematic_package

# 問題: 網路連線問題
# 解決方案: 更換源
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple/ package_name

# 問題: 權限錯誤
# 解決方案: 使用 --user 標記
pip install --user package_name
```

### R 套件問題
```bash
# 問題: Bioconductor 套件安裝失敗
# 解決方案: 先安裝 BiocManager
R -e "install.packages('BiocManager')"
R -e "BiocManager::install('problematic_package')"

# 問題: 編譯依賴缺失
# 解決方案: 安裝系統依賴
sudo apt-get update
sudo apt-get install r-base-dev
```

### 記憶體不足
```bash
# 問題: 大型套件安裝時記憶體不足
# 解決方案 1: 增加 swap 空間
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# 解決方案 2: 分批安裝
pip install large_package1
pip install large_package2  # 分開安裝

# 解決方案 3: 使用更輕量的替代品
pip install package_name_lite
```

## 📚 進階技巧

### 套件版本管理
```bash
# 鎖定套件版本
pip install "pandas==2.0.3"
echo "pandas==2.0.3" >> requirements.txt

# 版本範圍指定
pip install "numpy>=1.24.0,<1.26.0"

# 檢查可用版本
pip index versions pandas
```

### 開發模式安裝
```bash
# 可編輯安裝 (開發模式)
pip install -e /path/to/local/package

# 從 Git 倉庫安裝開發版
pip install git+https://github.com/user/repo.git@dev_branch
```

### 自動化腳本
```bash
#!/bin/bash
# auto-install.sh - 自動化套件安裝腳本

PACKAGES_FILE="my_packages.txt"

if [ -f "$PACKAGES_FILE" ]; then
    while read -r package; do
        echo "安裝套件: $package"
        bash .devcontainer/scripts/utils/package-manager.sh add python "$package"
    done < "$PACKAGES_FILE"
else
    echo "找不到套件清單檔案: $PACKAGES_FILE"
fi
```

## 🎯 最佳實踐

### 1. 環境管理
- **定期備份**: 重要變更前先備份
- **版本鎖定**: 生產環境使用固定版本
- **文檔記錄**: 記錄安裝的套件和原因

### 2. 套件選擇
- **官方優先**: 優先使用官方維護的套件
- **穩定版本**: 避免使用 alpha/beta 版本
- **相容性檢查**: 安裝前檢查與現有套件的相容性

### 3. 效能考量
- **批次操作**: 批次安裝多個套件
- **快取利用**: 妥善利用套件管理器快取
- **資源監控**: 監控安裝過程的資源使用

---

**透過強大的套件管理系統，讓您的分析環境始終保持最新且符合需求**  
*動態調整，無需重建，效率倍增！*