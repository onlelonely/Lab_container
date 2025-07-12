# 快速開始指南

**目標**: 5分鐘內啟動並運行您的第一個分析環境

## 🚀 5分鐘快速上手

### 步驟 1: 克隆專案
```bash
git clone https://github.com/onlelonely/Lab_container.git
cd Lab_container
```

### 步驟 2: 選擇分析環境
```bash
# 執行環境設置腳本
bash .devcontainer/scripts/manage/setup-environment.sh

# 互動式選擇菜單將出現:
# 1) minimal        - 基本 Python/R 環境
# 2) datascience    - 資料科學標準環境 
# 3) bioinformatics - NGS 生物資訊環境
# 4) ml             - 機器學習環境
# 5) statistics     - 統計分析環境
# 6) full           - 完整功能環境

# 建議新手選擇: 2 (datascience)
```

### 步驟 3: 開啟開發容器
```bash
# 在 VS Code 中開啟專案
code .

# VS Code 將自動檢測到 .devcontainer 設定
# 點擊右下角出現的 "Reopen in Container" 按鈕
```

### 步驟 4: 驗證環境
容器啟動後，在終端機中執行：
```bash
# 檢查 Python 環境
python --version
pip list

# 檢查 R 環境  
R --version
R -e "installed.packages()[,c('Package','Version')]"

# 檢查整體環境狀態
bash .devcontainer/scripts/manage/devcontainer-status.sh
```

## 🎯 常用環境快速設置

### 📊 資料科學環境 (推薦新手)
```bash
# 1. 設置環境
bash .devcontainer/scripts/manage/setup-environment.sh datascience

# 2. 驗證套件
python -c "import pandas, numpy, matplotlib; print('資料科學環境就緒!')"
R -e "library(tidyverse); cat('R 環境就緒!\n')"
```

**包含套件**:
- **Python**: pandas, numpy, matplotlib, seaborn, plotly, scikit-learn
- **R**: tidyverse, ggplot2, dplyr, readr
- **Jupyter**: JupyterLab 完整環境

### 🧬 生物資訊環境
```bash
# 1. 設置環境
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics

# 2. 驗證生物資訊工具
fastqc --version
samtools --version
```

**包含工具**:
- **輕量級**: FastQC, SAMtools, BCFtools, BEDtools
- **Python**: biopython, pysam, pyvcf
- **R**: Biostrings, GenomicRanges (透過 Bioconductor)

### 🤖 機器學習環境
```bash
# 1. 設置環境  
bash .devcontainer/scripts/manage/setup-environment.sh ml

# 2. 驗證 ML 套件
python -c "import sklearn, tensorflow; print('ML 環境就緒!')"
```

**包含套件**:
- **Python**: scikit-learn, tensorflow, torch, transformers
- **R**: caret, randomForest, nnet

## 🛠 環境管理快速命令

### 檢查環境狀態
```bash
# 快速狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh quick

# 完整狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh
```

### 動態套件管理
```bash
# 新增 Python 套件
bash .devcontainer/scripts/utils/package-manager.sh add python requests

# 新增 R 套件  
bash .devcontainer/scripts/utils/package-manager.sh add r forecast

# 列出已安裝套件
bash .devcontainer/scripts/utils/package-manager.sh list
```

### 環境備份
```bash
# 備份當前環境
bash .devcontainer/scripts/utils/package-manager.sh backup

# 備份將儲存在 backups/YYYYMMDD_HHMMSS/ 目錄
```

## 🧪 快速測試範例

### Python 資料分析範例
```python
# 在 Jupyter 或 Python 中執行
import pandas as pd
import matplotlib.pyplot as plt

# 建立測試資料
data = pd.DataFrame({
    'x': range(10),
    'y': [i**2 for i in range(10)]
})

# 繪圖
plt.plot(data['x'], data['y'])
plt.title('快速測試圖表')
plt.show()
```

### R 資料分析範例
```r
# 在 R 控制台或 RStudio 中執行
library(ggplot2)

# 建立測試資料
data <- data.frame(
  x = 1:10,
  y = (1:10)^2
)

# 繪圖
ggplot(data, aes(x, y)) + 
  geom_line() + 
  ggtitle("快速測試圖表")
```

### 生物資訊工具測試
```bash
# 測試 FastQC (需要有 FASTQ 檔案)
echo ">test_seq" > test.fq
echo "ATCGATCGATCGATCG" >> test.fq
fastqc test.fq

# 測試 SAMtools
echo "測試 SAMtools 安裝..."
samtools --help | head -5
```

## 🚨 常見問題快速解決

### 問題 1: 容器啟動失敗
```bash
# 解決方案: 重新生成設定檔
bash .devcontainer/scripts/manage/setup-environment.sh
```

### 問題 2: 套件缺失
```bash
# 解決方案: 重新安裝套件
bash .devcontainer/scripts/utils/package-manager.sh add python <套件名稱>
```

### 問題 3: 記憶體不足
```bash
# 解決方案: 檢查資源使用
bash .devcontainer/scripts/manage/devcontainer-status.sh
# 考慮使用 minimal 環境
```

## 📚 下一步學習

完成快速開始後，建議閱讀：

1. **[環境設置指南](./environment-setup.md)** - 了解詳細的環境配置選項
2. **[Profile 系統說明](./profiles.md)** - 深入了解環境配置檔案系統  
3. **[套件管理指南](./package-management.md)** - 學習進階套件管理技巧

## 🎉 成功標誌

如果您能夠：
- ✅ 成功啟動開發容器
- ✅ 執行 Python/R 程式碼
- ✅ 安裝新的套件
- ✅ 檢查環境狀態

**恭喜！您已經成功設置了 Lab Container 環境！**

---

**疑問或需要協助？**
- 查看 [常見問題解答](./faq.md)
- 提交 [GitHub Issue](https://github.com/onlelonely/Lab_container/issues)
- 參與 [社群討論](https://github.com/onlelonely/Lab_container/discussions)