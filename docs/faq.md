# 常見問題解答 (FAQ)

**文檔版本**: v2.0  
**適用系統**: Lab Container v2.0+  
**目標讀者**: 所有使用者

## 🚀 入門問題

### Q1: Lab Container 是什麼？
**A**: Lab Container 是一個基於 Docker 的科學計算環境，專門為資料科學、生物資訊學、機器學習等領域設計。它提供預配置的開發環境，讓您能夠快速開始分析工作，無需複雜的環境設置。

### Q2: 我需要什麼系統需求？
**A**: 
- **作業系統**: Windows 10/11, macOS 10.15+, 或 Linux (Ubuntu 18.04+)
- **記憶體**: 最少 4GB RAM (推薦 8GB+)
- **磁碟空間**: 最少 10GB 可用空間
- **軟體**:
  - Docker Desktop (Windows/macOS) 或 Docker Engine (Linux)
  - VS Code + Dev Containers Extension
  - Git

### Q3: 如何選擇適合的 Profile？
**A**: 根據您的主要工作領域選擇：
- **新手學習**: `minimal` - 輕量化環境
- **資料分析**: `datascience` - Python + R 資料科學工具
- **基因分析**: `bioinformatics` - NGS 資料分析工具
- **機器學習**: `ml` - 深度學習框架
- **統計分析**: `statistics` - R 統計套件
- **多領域**: `full` - 完整功能

### Q4: 第一次啟動需要多長時間？
**A**:
- **minimal**: 5-10 分鐘
- **datascience**: 10-15 分鐘  
- **bioinformatics**: 15-20 分鐘
- **ml**: 20-30 分鐘
- **statistics**: 10-15 分鐘
- **full**: 25-35 分鐘

*實際時間取決於網路速度和系統效能*

## 🔧 使用問題

### Q5: 如何安裝新的 Python/R 套件？
**A**: 使用內建的套件管理器：
```bash
# Python 套件
bash .devcontainer/scripts/utils/package-manager.sh add python pandas

# R 套件
bash .devcontainer/scripts/utils/package-manager.sh add r ggplot2

# Conda 套件
bash .devcontainer/scripts/utils/package-manager.sh add conda numpy
```

### Q6: 我可以在運行時切換 Profile 嗎？
**A**: 是的，但需要重建容器：
```bash
bash .devcontainer/scripts/manage/setup-environment.sh <new_profile>
```
*注意：這會重置環境，建議先備份自訂的套件*

### Q7: 如何備份我的環境配置？
**A**: 使用自動備份功能：
```bash
# 建立備份
bash .devcontainer/scripts/utils/package-manager.sh backup

# 還原備份
bash .devcontainer/scripts/utils/package-manager.sh restore latest
```

### Q8: 如何使用重量級生物資訊工具 (如 GATK)？
**A**: 透過工具管理器：
```bash
# 檢查工具可用性
bash .devcontainer/scripts/manage/tool-manager.sh check gatk

# 執行 GATK 命令
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller --help
```

## 🐛 問題解決

### Q9: 容器啟動失敗怎麼辦？
**A**: 按以下順序檢查：
1. 確認 Docker 服務正在運行：`docker version`
2. 檢查磁碟空間：`df -h`
3. 查看錯誤日誌：`docker-compose logs`
4. 清理 Docker 快取：`docker system prune`
5. 重新建構：在 VS Code 中選擇 "Rebuild Container"

### Q10: 套件安裝失敗如何解決？
**A**: 常見解決方法：
```bash
# Python 套件編譯錯誤 → 使用 Conda
bash .devcontainer/scripts/utils/package-manager.sh add conda problem_package

# R 套件缺少依賴 → 安裝系統依賴
sudo apt-get install r-base-dev libxml2-dev

# 網路超時 → 使用國內源
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple/ package_name
```

### Q11: GPU 不被識別怎麼辦？
**A**: 檢查 GPU 支援：
```bash
# 檢查 NVIDIA 驅動
nvidia-smi

# 檢查 Docker GPU 支援
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi

# 重新安裝 NVIDIA Container Toolkit
sudo apt-get install nvidia-container-toolkit
sudo systemctl restart docker
```

### Q12: 記憶體不足錯誤如何處理？
**A**:
1. **增加虛擬記憶體**：
```bash
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

2. **切換到輕量 Profile**：
```bash
bash .devcontainer/scripts/manage/setup-environment.sh minimal
```

3. **分批安裝大型套件**

## 🔄 工作流程

### Q13: 如何進行 NGS 資料分析？
**A**: 選擇 `bioinformatics` Profile，然後：
```bash
# 1. 品質控制
fastqc data/*.fastq.gz -o results/qc/

# 2. 序列比對
bwa mem reference.fa reads.fastq.gz | samtools sort -o aligned.bam

# 3. 變異呼叫
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I aligned.bam -O variants.vcf -R reference.fa
```

### Q14: 如何開始機器學習專案？
**A**: 選擇 `ml` Profile，然後：
```bash
# 1. 啟動 Jupyter Lab
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser

# 2. 檢查 GPU 可用性
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# 3. 開始訓練模型
```

### Q15: 如何生成統計報告？
**A**: 選擇 `statistics` Profile，使用 R Markdown：
```r
# 在 R 中建立新報告
rmarkdown::render("analysis.Rmd", output_format = "html_document")
```

## 🛠 自訂與擴展

### Q16: 如何建立自訂 Profile？
**A**:
```bash
# 1. 複製現有 Profile
cp .devcontainer/profiles/datascience.env .devcontainer/profiles/myprofile.env

# 2. 編輯配置
vim .devcontainer/profiles/myprofile.env

# 3. 建立套件清單
touch .devcontainer/configs/python/requirements-myprofile.txt

# 4. 測試新 Profile
bash .devcontainer/scripts/manage/setup-environment.sh myprofile
```

### Q17: 如何新增自訂工具？
**A**: 
**輕量級工具** (整合到主容器)：
```bash
# 更新 Conda 環境檔案
vim .devcontainer/configs/conda/environment-<profile>.yml
```

**重量級工具** (獨立容器)：
```bash
# 更新 docker-compose.tools.yml
vim .devcontainer/docker-compose.tools.yml
```

### Q18: 如何設定自動啟動的服務？
**A**: 編輯 `devcontainer.json`：
```json
{
  "postCreateCommand": "bash .devcontainer/scripts/setup/install-additional-tools.sh",
  "postStartCommand": "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser &"
}
```

## 📊 效能最佳化

### Q19: 如何加速容器啟動？
**A**:
1. **使用本地快取**：避免重複下載
2. **選擇合適的 Profile**：不要使用過於複雜的環境
3. **清理無用資源**：`docker system prune`
4. **使用 SSD 硬碟**：顯著提升 I/O 效能

### Q20: 如何優化大資料集處理？
**A**:
```bash
# 1. 使用記憶體映射
# Python 中使用 mmap 或 pandas chunking
pandas.read_csv('large_file.csv', chunksize=10000)

# 2. 使用壓縮格式
# 偏好 parquet 格式而非 CSV

# 3. 增加暫存空間
export TMPDIR=/path/to/large/tmp

# 4. 調整 Docker 資源限制
# 編輯 docker-compose.yml 增加記憶體限制
```

## 🔒 安全性問題

### Q21: 容器內的資料是否安全？
**A**: 
- 資料透過 Docker Volume 掛載，儲存在主機系統
- 容器刪除不會影響掛載的資料
- 建議定期備份重要資料到外部儲存

### Q22: 如何保護敏感資料？
**A**:
```bash
# 1. 使用環境變數儲存敏感資訊
echo "API_KEY=your_secret_key" >> .env

# 2. 避免將敏感檔案加入版本控制
echo "*.key" >> .gitignore
echo ".env" >> .gitignore

# 3. 使用加密儲存
gpg --symmetric sensitive_data.txt
```

## 🌐 協作與分享

### Q23: 如何與團隊分享環境配置？
**A**:
```bash
# 1. 匯出環境配置
bash .devcontainer/scripts/utils/package-manager.sh backup

# 2. 版本控制 Profile 設定
git add .devcontainer/profiles/
git commit -m "Add team profile configuration"

# 3. 分享配置檔案
# 團隊成員可以使用相同的 Profile 重建環境
```

### Q24: 如何在不同機器間同步環境？
**A**:
```bash
# 機器 A - 匯出環境
bash .devcontainer/scripts/utils/package-manager.sh backup

# 機器 B - 匯入環境  
bash .devcontainer/scripts/utils/package-manager.sh restore path/to/backup
```

## 📚 學習資源

### Q25: 在哪裡可以找到更多學習資料？
**A**:
- **官方文檔**: `docs/` 目錄下的完整文檔
- **範例腳本**: `.devcontainer/examples/` 目錄
- **工作流程範例**: `docs/workflows.md`
- **故障排除**: `docs/troubleshooting.md`

### Q26: 如何貢獻代碼或報告問題？
**A**:
- **報告 Bug**: 在 GitHub Issues 中建立問題報告
- **功能請求**: 提交功能需求說明
- **代碼貢獻**: Fork → 修改 → Pull Request
- **文檔改進**: 直接提交文檔更新的 PR

## 🆕 版本更新

### Q27: 如何更新到新版本？
**A**:
```bash
# 1. 備份當前環境
bash .devcontainer/scripts/utils/package-manager.sh backup

# 2. 拉取最新代碼
git pull origin main

# 3. 重建容器
# 在 VS Code 中選擇 "Rebuild Container"

# 4. 如有問題，還原備份
bash .devcontainer/scripts/utils/package-manager.sh restore latest
```

### Q28: 新版本是否向後相容？
**A**: 
- **配置檔案**: 向後相容，舊 Profile 仍可使用
- **API 介面**: 主要命令保持穩定
- **資料格式**: 備份格式向後相容
- **重大變更**: 會在發行說明中明確標註

---

**如果您的問題未在此列出，請查閱詳細文檔或建立 GitHub Issue**  
*我們會持續更新此 FAQ 以涵蓋更多使用場景*