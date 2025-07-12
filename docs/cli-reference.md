# 命令列參考手冊

**文檔版本**: v2.0  
**適用系統**: Lab Container v2.0+  
**目標讀者**: 所有使用者

## 📚 命令列工具概覽

Lab Container 提供完整的命令列工具集，讓您能夠透過腳本和終端高效管理整個環境。

### 工具分類
```
命令列工具結構:
├── 🛠️ 環境管理工具
│   ├── setup-environment.sh      # 環境設置與配置
│   └── devcontainer-status.sh    # 系統狀態檢查
├── 📦 套件管理工具  
│   └── package-manager.sh         # 套件安裝、移除、備份
├── 🧬 工具管理器
│   └── tool-manager.sh            # 生物資訊工具管理
└── 🔧 輔助工具
    ├── backup-restore.sh          # 備份與還原
    └── cleanup.sh                 # 清理與最佳化
```

## 🛠 環境管理命令

### setup-environment.sh
**功能**: 環境設置與 Profile 管理

#### 基本語法
```bash
bash .devcontainer/scripts/manage/setup-environment.sh [PROFILE] [OPTIONS]
```

#### 參數說明
```bash
PROFILE                    # 選擇的環境配置檔案
  minimal                  # 最小化環境
  datascience             # 資料科學環境  
  bioinformatics          # 生物資訊環境
  ml                      # 機器學習環境
  statistics              # 統計分析環境
  full                    # 完整功能環境

OPTIONS
  --force                 # 強制重新生成配置
  --dry-run              # 預覽配置但不執行
  --verbose              # 詳細輸出模式
  --no-backup            # 跳過自動備份
  --help                 # 顯示說明
```

#### 使用範例
```bash
# 互動式選擇 Profile
bash .devcontainer/scripts/manage/setup-environment.sh

# 直接指定 Profile
bash .devcontainer/scripts/manage/setup-environment.sh datascience

# 強制重新配置
bash .devcontainer/scripts/manage/setup-environment.sh full --force

# 預覽配置變更
bash .devcontainer/scripts/manage/setup-environment.sh ml --dry-run

# 詳細模式設置
bash .devcontainer/scripts/manage/setup-environment.sh bioinformatics --verbose
```

### devcontainer-status.sh
**功能**: 系統狀態檢查與診斷

#### 基本語法
```bash
bash .devcontainer/scripts/manage/devcontainer-status.sh [OPTIONS]
```

#### 參數說明
```bash
OPTIONS
  --verbose              # 詳細狀態報告
  --component COMP       # 檢查特定組件
  --format FORMAT        # 輸出格式 (text|json|html)
  --output FILE          # 輸出到檔案
  --quiet               # 僅顯示錯誤
  --help                # 顯示說明

COMPONENTS
  docker                # Docker 服務狀態
  containers            # 容器狀態
  tools                 # 工具可用性
  packages              # 套件狀態
  network               # 網路連線
  storage               # 儲存空間
```

#### 使用範例
```bash
# 基本狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh

# 詳細狀態報告
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose

# 檢查特定組件
bash .devcontainer/scripts/manage/devcontainer-status.sh --component docker
bash .devcontainer/scripts/manage/devcontainer-status.sh --component tools

# JSON 格式輸出
bash .devcontainer/scripts/manage/devcontainer-status.sh --format json

# 輸出到檔案
bash .devcontainer/scripts/manage/devcontainer-status.sh --output status_report.html --format html
```

## 📦 套件管理命令

### package-manager.sh
**功能**: 動態套件管理系統

#### 基本語法
```bash
bash .devcontainer/scripts/utils/package-manager.sh <ACTION> [TYPE] [PACKAGE] [OPTIONS]
```

#### 動作 (ACTION)
```bash
add                       # 安裝套件
remove                    # 移除套件
list                      # 列出已安裝套件
search                    # 搜尋套件
update                    # 更新套件
backup                    # 備份環境
restore                   # 還原環境
check                     # 檢查套件狀態
help                      # 顯示說明
```

#### 套件類型 (TYPE)
```bash
python                    # Python 套件 (pip)
r                        # R 套件 (CRAN/Bioconductor)
conda                    # Conda 套件
system                   # 系統套件 (apt)
```

#### 選項 (OPTIONS)
```bash
--version VERSION        # 指定版本
--channel CHANNEL        # 指定 Conda channel
--source SOURCE          # 指定套件源
--force                  # 強制執行
--dry-run               # 預覽但不執行
--verbose               # 詳細輸出
--quiet                 # 靜默模式
```

#### 使用範例

##### 安裝套件
```bash
# Python 套件
bash .devcontainer/scripts/utils/package-manager.sh add python pandas
bash .devcontainer/scripts/utils/package-manager.sh add python "tensorflow>=2.13.0"
bash .devcontainer/scripts/utils/package-manager.sh add python "torch==2.0.1" --force

# R 套件
bash .devcontainer/scripts/utils/package-manager.sh add r ggplot2
bash .devcontainer/scripts/utils/package-manager.sh add r "dplyr tidyr readr"

# Conda 套件
bash .devcontainer/scripts/utils/package-manager.sh add conda numpy
bash .devcontainer/scripts/utils/package-manager.sh add conda "bioconda::samtools" --channel bioconda

# 系統套件
bash .devcontainer/scripts/utils/package-manager.sh add system curl wget git
```

##### 移除套件
```bash
# 移除 Python 套件
bash .devcontainer/scripts/utils/package-manager.sh remove python pandas

# 移除多個套件
bash .devcontainer/scripts/utils/package-manager.sh remove python "numpy matplotlib"

# 強制移除
bash .devcontainer/scripts/utils/package-manager.sh remove python tensorflow --force
```

##### 列出和搜尋
```bash
# 列出所有套件
bash .devcontainer/scripts/utils/package-manager.sh list

# 列出特定類型套件
bash .devcontainer/scripts/utils/package-manager.sh list python
bash .devcontainer/scripts/utils/package-manager.sh list r

# 搜尋套件
bash .devcontainer/scripts/utils/package-manager.sh search python tensorflow
bash .devcontainer/scripts/utils/package-manager.sh search conda samtools
```

##### 備份與還原
```bash
# 建立環境備份
bash .devcontainer/scripts/utils/package-manager.sh backup

# 備份到指定位置
bash .devcontainer/scripts/utils/package-manager.sh backup --output /path/to/backup

# 還原最新備份
bash .devcontainer/scripts/utils/package-manager.sh restore latest

# 還原特定備份
bash .devcontainer/scripts/utils/package-manager.sh restore backups/20250712_143000

# 還原時跳過確認
bash .devcontainer/scripts/utils/package-manager.sh restore latest --force
```

## 🧬 工具管理命令

### tool-manager.sh
**功能**: 生物資訊工具容器管理

#### 基本語法
```bash
bash .devcontainer/scripts/manage/tool-manager.sh <ACTION> [TOOL] [ARGS...]
```

#### 動作 (ACTION)
```bash
# 工具檢查
check                     # 檢查工具可用性
health                    # 健康狀態檢查

# 容器管理
manage                    # 容器生命週期管理
  start                   # 啟動容器
  stop                    # 停止容器
  restart                 # 重啟容器
  status                  # 容器狀態
  logs                    # 查看日誌
  clean                   # 清理容器

# 工具執行
gatk                     # 執行 GATK 命令
colabfold                # 執行 ColabFold 命令
autodock                 # 執行 AutoDock Vina 命令

# 輔助功能
help                     # 顯示說明
version                  # 顯示版本資訊
```

#### 使用範例

##### 工具檢查
```bash
# 檢查所有工具
bash .devcontainer/scripts/manage/tool-manager.sh check all

# 檢查特定工具
bash .devcontainer/scripts/manage/tool-manager.sh check gatk
bash .devcontainer/scripts/manage/tool-manager.sh check colabfold

# 健康狀態檢查
bash .devcontainer/scripts/manage/tool-manager.sh health
```

##### 容器管理
```bash
# 啟動工具容器
bash .devcontainer/scripts/manage/tool-manager.sh manage start gatk
bash .devcontainer/scripts/manage/tool-manager.sh manage start colabfold

# 停止容器
bash .devcontainer/scripts/manage/tool-manager.sh manage stop gatk
bash .devcontainer/scripts/manage/tool-manager.sh manage stop all

# 檢查狀態
bash .devcontainer/scripts/manage/tool-manager.sh manage status

# 查看日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk -f  # 即時日誌

# 清理未使用容器
bash .devcontainer/scripts/manage/tool-manager.sh manage clean
```

##### 工具執行
```bash
# GATK 命令
bash .devcontainer/scripts/manage/tool-manager.sh gatk --version
bash .devcontainer/scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I input.bam -O variants.vcf -R reference.fa

# ColabFold 結構預測
bash .devcontainer/scripts/manage/tool-manager.sh colabfold \
    protein.fasta results/ --num-models 5

# AutoDock Vina 分子對接
bash .devcontainer/scripts/manage/tool-manager.sh autodock \
    --receptor receptor.pdbqt --ligand ligand.pdbqt --out result.pdbqt
```

## 🔧 輔助工具命令

### backup-restore.sh
**功能**: 進階備份與還原功能

#### 基本語法
```bash
bash .devcontainer/scripts/utils/backup-restore.sh <ACTION> [OPTIONS]
```

#### 動作與選項
```bash
# 備份動作
backup
  --full                  # 完整系統備份
  --config-only           # 僅備份配置
  --data-only            # 僅備份資料
  --compress             # 壓縮備份
  --encrypt              # 加密備份

# 還原動作  
restore
  --backup-path PATH     # 指定備份路徑
  --selective ITEMS      # 選擇性還原
  --verify               # 驗證還原完整性

# 清理動作
cleanup
  --old-backups          # 清理舊備份
  --temp-files           # 清理暫存檔案
  --logs                 # 清理日誌檔案
```

#### 使用範例
```bash
# 完整系統備份
bash .devcontainer/scripts/utils/backup-restore.sh backup --full --compress

# 僅備份配置檔案
bash .devcontainer/scripts/utils/backup-restore.sh backup --config-only

# 加密備份
bash .devcontainer/scripts/utils/backup-restore.sh backup --full --encrypt

# 選擇性還原
bash .devcontainer/scripts/utils/backup-restore.sh restore \
    --backup-path backups/20250712_143000 \
    --selective "configs,packages"

# 清理舊備份
bash .devcontainer/scripts/utils/backup-restore.sh cleanup --old-backups
```

### cleanup.sh
**功能**: 系統清理與最佳化

#### 基本語法
```bash
bash .devcontainer/scripts/utils/cleanup.sh [OPTIONS]
```

#### 選項說明
```bash
--docker                 # 清理 Docker 資源
--packages               # 清理套件快取
--logs                   # 清理日誌檔案
--temp                   # 清理暫存檔案
--all                    # 全面清理
--dry-run               # 預覽清理項目
--force                 # 強制清理 (不詢問)
--verbose               # 詳細輸出
```

#### 使用範例
```bash
# 基本清理
bash .devcontainer/scripts/utils/cleanup.sh

# Docker 資源清理
bash .devcontainer/scripts/utils/cleanup.sh --docker

# 全面清理
bash .devcontainer/scripts/utils/cleanup.sh --all

# 預覽清理項目
bash .devcontainer/scripts/utils/cleanup.sh --all --dry-run

# 強制清理 (不詢問確認)
bash .devcontainer/scripts/utils/cleanup.sh --all --force
```

## 🔄 工作流程命令組合

### 常用命令組合範例

#### 環境初始化
```bash
# 完整環境設置流程
bash .devcontainer/scripts/manage/setup-environment.sh datascience
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose
bash .devcontainer/scripts/utils/package-manager.sh list
```

#### 日常維護
```bash
# 每日維護流程
bash .devcontainer/scripts/manage/devcontainer-status.sh --quiet
bash .devcontainer/scripts/utils/cleanup.sh --temp --logs
bash .devcontainer/scripts/utils/package-manager.sh backup
```

#### 問題診斷
```bash
# 系統診斷流程
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose --output diagnosis.html
bash .devcontainer/scripts/manage/tool-manager.sh check all
docker logs lab-container-devcontainer-1 --tail 100
```

#### 環境遷移
```bash
# 環境匯出
bash .devcontainer/scripts/utils/package-manager.sh backup --output migration_backup
bash .devcontainer/scripts/utils/backup-restore.sh backup --full --compress

# 環境匯入
bash .devcontainer/scripts/utils/package-manager.sh restore migration_backup
bash .devcontainer/scripts/manage/setup-environment.sh <profile> --force
```

## 📋 全域選項說明

### 通用選項
大部分腳本都支援以下通用選項：

```bash
--help                   # 顯示命令說明
--version               # 顯示版本資訊
--verbose               # 詳細輸出模式
--quiet                 # 靜默模式 (僅錯誤)
--dry-run              # 預覽模式 (不執行實際操作)
--force                # 強制執行 (跳過確認)
--config CONFIG_FILE   # 指定配置檔案
--log-file LOG_FILE    # 指定日誌檔案
```

### 環境變數
```bash
# 可用的環境變數
export DEVCONTAINER_DEBUG=1           # 啟用除錯模式
export DEVCONTAINER_NO_BACKUP=1       # 跳過自動備份
export DEVCONTAINER_PARALLEL_JOBS=4   # 並行工作數
export DEVCONTAINER_TIMEOUT=300       # 命令超時時間 (秒)
export DEVCONTAINER_CACHE_DIR=/tmp    # 快取目錄
```

## 🎯 最佳實踐

### 1. 腳本使用建議
- **始終閱讀 --help**: 每個腳本都有詳細說明
- **使用 --dry-run**: 重要操作前先預覽
- **定期備份**: 重大變更前先備份
- **記錄日誌**: 重要操作時指定日誌檔案

### 2. 常見模式
```bash
# 安全的套件安裝模式
bash .devcontainer/scripts/utils/package-manager.sh backup
bash .devcontainer/scripts/utils/package-manager.sh add python package_name --dry-run
bash .devcontainer/scripts/utils/package-manager.sh add python package_name

# 系統維護模式
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose
bash .devcontainer/scripts/utils/cleanup.sh --dry-run
bash .devcontainer/scripts/utils/cleanup.sh
```

### 3. 故障排除
```bash
# 當命令失敗時的診斷步驟
bash .devcontainer/scripts/manage/devcontainer-status.sh --component <relevant_component>
docker logs lab-container-devcontainer-1 --tail 50
bash <script_name> --verbose --dry-run
```

---

**掌握這些命令列工具，讓您成為 Lab Container 的專家使用者**  
*功能強大，使用簡單，效率倍增！*