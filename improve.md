# Dev Container 彈性環境改進建議

**撰寫日期**: 2025-07-11  
**狀態**: 等待 CI/CD 完成後實施  
**目標**: 提升 dev container 的彈性，支援各種分析環境組合

## 📋 當前架構評估

### ✅ 已具備的優勢
1. **完整的重構基礎** - 已按照 `refactor.md` 完成系統性重構
2. **多階段構建** - 支援 core/extended 不同功能組合
3. **完整的語言支援** - Python 3.11, R 4.3, Jupyter Lab
4. **NGS 生物資訊工具** - Bioconductor, BioPython 等完整工具鏈
5. **模組化配置** - 分離的 requirements 和 conda 環境檔案
6. **穩定的錯誤處理** - 重試邏輯和驗證機制
7. **安全性** - 非根用戶配置，完整的權限管理

### 🎯 改進目標
基於您的需求「彈性組成各種分析環境」，當前系統已經很好，但可以在以下方面增強：

1. **環境選擇的簡化** - 讓使用者更容易選擇所需的環境組合
2. **動態套件管理** - 支援運行時新增/移除套件
3. **預設環境模板** - 提供常見的分析環境組合
4. **配置的可視化** - 更好的配置管理界面

## 🛠 具體改進建議

### 1. 環境設定檔(Profile)系統

#### 1.1 建立預設 Profile
在 `.devcontainer/profiles/` 目錄建立不同的環境設定檔：

```bash
profiles/
├── minimal.env           # 基本 Python/R 環境
├── datascience.env      # 資料科學標準環境
├── bioinformatics.env   # NGS 生物資訊環境
├── ml.env               # 機器學習環境
├── statistics.env       # 統計分析環境
└── full.env             # 完整功能環境
```

#### 1.2 Profile 內容範例
```bash
# profiles/datascience.env
PROFILE_NAME="Data Science"
PROFILE_DESCRIPTION="標準資料科學環境，包含 pandas, numpy, matplotlib, seaborn"

# Python 套件
PYTHON_PACKAGES="pandas numpy matplotlib seaborn plotly scikit-learn"
PYTHON_REQUIREMENTS_FILE="requirements-datascience.txt"

# R 套件
R_PACKAGES="tidyverse ggplot2 dplyr readr"
R_SCRIPT_FILE="datascience-packages.R"

# Conda 環境
CONDA_ENV_FILE="environment-datascience.yml"

# 啟用的功能
ENABLE_JUPYTER=true
ENABLE_RSTUDIO=false
ENABLE_VSCODE_R=true
ENABLE_VSCODE_PYTHON=true
```

### 2. 智能環境選擇器

#### 2.1 互動式環境選擇
建立 `setup-environment.sh` 腳本：

```bash
#!/bin/bash
# setup-environment.sh

echo "🔧 Dev Container 環境設定"
echo "請選擇您需要的分析環境："
echo ""
echo "1) minimal        - 基本 Python/R 環境"
echo "2) datascience    - 資料科學標準環境 (pandas, tidyverse)"
echo "3) bioinformatics - NGS 生物資訊環境 (Bioconductor, BioPython)"
echo "4) ml             - 機器學習環境 (scikit-learn, tensorflow)"
echo "5) statistics     - 統計分析環境 (R focus)"
echo "6) full           - 完整功能環境"
echo "7) custom         - 自訂環境"
echo ""
read -p "請輸入選項 (1-7): " choice

case $choice in
    1) load_profile "minimal" ;;
    2) load_profile "datascience" ;;
    3) load_profile "bioinformatics" ;;
    4) load_profile "ml" ;;
    5) load_profile "statistics" ;;
    6) load_profile "full" ;;
    7) setup_custom_environment ;;
    *) echo "無效選項，使用預設環境" && load_profile "datascience" ;;
esac
```

#### 2.2 Profile 載入系統
```bash
load_profile() {
    local profile_name="$1"
    local profile_file=".devcontainer/profiles/${profile_name}.env"
    
    if [ -f "$profile_file" ]; then
        source "$profile_file"
        echo "✅ 載入環境設定: $PROFILE_NAME"
        generate_devcontainer_config
        generate_docker_compose_config
    else
        echo "❌ 找不到環境設定檔: $profile_file"
        exit 1
    fi
}
```

### 3. 動態套件管理系統

#### 3.1 運行時套件管理
擴展 `config-manager.sh` 支援動態套件管理：

```bash
# 新增功能
devcontainer-add-package() {
    local package_type="$1"  # python, r, conda
    local package_name="$2"
    
    case $package_type in
        python)
            pip install "$package_name"
            echo "$package_name" >> .devcontainer/configs/python/requirements-runtime.txt
            ;;
        r)
            R -e "install.packages('$package_name')"
            echo "install.packages('$package_name')" >> .devcontainer/configs/r/runtime-packages.R
            ;;
        conda)
            conda install -y "$package_name"
            echo "  - $package_name" >> .devcontainer/configs/conda/runtime-packages.yml
            ;;
    esac
}

# 使用範例
devcontainer-add-package python "requests"
devcontainer-add-package r "forecast"
devcontainer-add-package conda "nodejs"
```

#### 3.2 套件狀態管理
```bash
# 查看已安裝套件
devcontainer-list-packages() {
    echo "=== Python 套件 ==="
    pip list
    echo ""
    echo "=== R 套件 ==="
    R -e "installed.packages()[,c('Package','Version')]"
    echo ""
    echo "=== Conda 套件 ==="
    conda list
}

# 備份當前環境
devcontainer-export-environment() {
    local backup_dir="backups/$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$backup_dir"
    
    pip freeze > "$backup_dir/requirements.txt"
    R -e "installed.packages()[,c('Package','Version')] |> write.csv('$backup_dir/r-packages.csv')"
    conda env export > "$backup_dir/environment.yml"
    
    echo "環境已備份至: $backup_dir"
}
```

### 4. 配置範本系統

#### 4.1 動態生成 devcontainer.json
```bash
generate_devcontainer_config() {
    local template_file=".devcontainer/templates/devcontainer.json.template"
    local output_file=".devcontainer/devcontainer.json"
    
    # 替換範本中的變數
    envsubst < "$template_file" > "$output_file"
    
    echo "✅ 已生成 devcontainer.json"
}
```

#### 4.2 範本檔案範例
```json
{
  "name": "${PROFILE_NAME}",
  "dockerComposeFile": "docker-compose.yml",
  "service": "devcontainer",
  "workspaceFolder": "/workspace",
  "shutdownAction": "stopCompose",
  "customizations": {
    "vscode": {
      "extensions": [
        "ms-python.python",
        "ms-python.black-formatter",
        "ms-toolsai.jupyter"
      ]
    }
  },
  "remoteEnv": {
    "PROFILE": "${PROFILE_NAME}",
    "PYTHON_VERSION": "${PYTHON_VERSION}",
    "R_VERSION": "${R_VERSION}"
  }
}
```

### 5. 使用者體驗改進

#### 5.1 環境狀態儀表板
建立 `devcontainer-status.sh`:

```bash
#!/bin/bash
# devcontainer-status.sh

echo "📊 Dev Container 環境狀態"
echo "=========================="
echo ""
echo "🔧 當前環境: $(cat .devcontainer/.current-profile 2>/dev/null || echo '未設定')"
echo "🐍 Python: $(python --version 2>/dev/null || echo '未安裝')"
echo "📊 R: $(R --version | head -1 2>/dev/null || echo '未安裝')"
echo "📚 Jupyter: $(jupyter --version 2>/dev/null || echo '未安裝')"
echo ""
echo "💾 磁碟使用:"
df -h /workspace 2>/dev/null || echo "無法取得磁碟資訊"
echo ""
echo "📦 套件統計:"
echo "   Python: $(pip list 2>/dev/null | wc -l || echo '0') 個套件"
echo "   R: $(R -e 'length(installed.packages()[,1])' 2>/dev/null | tail -1 || echo '0 個套件')"
echo ""
echo "🕒 最後更新: $(date)"
```

#### 5.2 快速啟動腳本
```bash
#!/bin/bash
# quick-start.sh

echo "🚀 Dev Container 快速啟動"
echo ""

# 檢查是否已有設定
if [ -f ".devcontainer/.current-profile" ]; then
    current_profile=$(cat .devcontainer/.current-profile)
    echo "發現現有環境: $current_profile"
    read -p "是否要重新設定環境? (y/N): " reset_env
    
    if [[ ! "$reset_env" =~ ^[Yy]$ ]]; then
        echo "使用現有環境啟動..."
        code .
        exit 0
    fi
fi

# 執行環境設定
bash .devcontainer/scripts/manage/setup-environment.sh

# 啟動 VS Code
echo "正在啟動 VS Code..."
code .
```

## 🗂 新增檔案結構

```
.devcontainer/
├── profiles/                    # 環境設定檔
│   ├── minimal.env
│   ├── datascience.env
│   ├── bioinformatics.env
│   ├── ml.env
│   ├── statistics.env
│   └── full.env
├── templates/                   # 配置範本
│   ├── devcontainer.json.template
│   └── docker-compose.yml.template
├── scripts/
│   ├── manage/
│   │   ├── setup-environment.sh      # 環境選擇器
│   │   ├── devcontainer-status.sh    # 狀態儀表板
│   │   └── quick-start.sh           # 快速啟動
│   └── utils/
│       ├── package-manager.sh       # 動態套件管理
│       └── profile-loader.sh        # Profile 載入器
├── configs/
│   ├── python/
│   │   └── requirements-runtime.txt  # 運行時安裝的套件
│   ├── r/
│   │   └── runtime-packages.R       # 運行時安裝的 R 套件
│   └── conda/
│       └── runtime-packages.yml     # 運行時安裝的 conda 套件
└── backups/                         # 環境備份
```

## 📅 實施計劃

### 階段一：基礎 Profile 系統 (CI/CD 完成後)
1. 建立 profiles 目錄和基本設定檔
2. 實作 profile 載入系統
3. 建立環境選擇器腳本
4. 測試基本的環境切換功能

### 階段二：動態套件管理
1. 擴展 config-manager.sh 支援動態套件管理
2. 建立套件狀態追蹤系統
3. 實作環境備份和還原功能
4. 測試套件的新增和移除

### 階段三：使用者體驗優化
1. 建立狀態儀表板
2. 實作快速啟動腳本
3. 改進錯誤訊息和使用者指引
4. 建立完整的使用說明

### 階段四：進階功能
1. 實作配置範本系統
2. 支援自訂環境建立
3. 建立環境分享和匯入功能
4. 整合 CI/CD 自動化測試

## 🔍 技術考量

### 1. 向後相容性
- 保持現有的環境設定不受影響
- 新功能為可選的增強功能
- 提供遷移路徑給現有使用者

### 2. 效能影響
- Profile 載入應該要快速
- 套件管理不應影響現有環境
- 備份和還原功能要高效

### 3. 維護性
- 保持模組化設計
- 清楚的設定檔案結構
- 完整的錯誤處理和日誌

## 🎯 預期效益

1. **提升彈性** - 使用者可以輕鬆選擇和切換不同的分析環境
2. **簡化管理** - 透過 Profile 系統減少配置複雜度
3. **動態調整** - 支援運行時新增套件，無需重建容器
4. **更好的使用者體驗** - 互動式設定和狀態儀表板
5. **環境一致性** - 標準化的環境配置和管理

## 📝 後續行動

1. **等待 CI/CD 完成** - 確保當前重構穩定
2. **評估實施順序** - 根據實際需求調整階段計劃
3. **建立測試計劃** - 確保新功能不影響現有功能
4. **撰寫詳細文檔** - 為每個新功能建立使用指南

---

**備註**: 此改進計劃建立在現有優秀架構的基礎上，目標是增強而非重構。所有建議都保持向後相容性，並且可以分階段實施。

---

# 🧬 專業生物資訊工具整合策略

**新增日期**: 2025-07-11  
**目標**: 整合專業生物資訊工具 (FastQC, SAMTools, GATK, LocalColabFold, AutoDock Vina)

## 🔍 工具特性分析與架構決策

### 輕量級工具 (整合到 Dev Container)
適合納入現有 dev container 的工具：

| 工具 | 特性 | 記憶體需求 | 儲存需求 | 使用頻率 |
|------|------|------------|----------|----------|
| **FastQC** | Java 小工具，快速 QC 分析 | < 1GB | < 100MB | 高 |
| **SAMTools** | C 工具，基本 BAM/SAM 操作 | < 2GB | < 500MB | 高 |
| **BCFtools** | VCF 檔案處理 | < 1GB | < 200MB | 高 |
| **BEDtools** | 基因組區間操作 | < 1GB | < 100MB | 高 |

### 重量級工具 (專用容器)
建議使用專用容器的工具：

| 工具 | 特性 | 記憶體需求 | 儲存需求 | GPU 需求 |
|------|------|------------|----------|----------|
| **GATK** | 大型變異呼叫工具包 | 8-32GB | 1-2GB | 無 |
| **LocalColabFold** | 蛋白質結構預測 | 16-64GB | 50-100GB | 建議 |
| **AutoDock Vina** | 分子對接工具 | 4-8GB | 500MB | 無 |

## 🏗 推薦架構策略：混合容器系統

### 架構圖
```
📦 Main Dev Container (現有)
├── 🐍 Python/R 分析環境
├── 📊 Jupyter Lab/RStudio
├── 🧬 輕量級生物資訊工具
│   ├── FastQC
│   ├── SAMTools
│   ├── BCFtools
│   └── BEDtools
└── 📝 資料分析腳本

📦 專用工具容器群
├── 🔬 GATK Container (變異呼叫)
├── 🧬 ColabFold Container (結構預測)
└── 💊 AutoDock Container (分子對接)
```

## 🛠 具體實施方案

### 1. 輕量級工具整合 (Phase 1)

#### 1.1 更新 bioinformatics profile
```bash
# profiles/bioinformatics.env
PROFILE_NAME="Bioinformatics"
PROFILE_DESCRIPTION="NGS 生物資訊分析環境，包含常用工具"

# 輕量級工具
BIOINFORMATICS_TOOLS="fastqc samtools bcftools bedtools"

# Python 生物資訊套件
PYTHON_PACKAGES="biopython pysam pyvcf"
PYTHON_REQUIREMENTS_FILE="requirements-bioinformatics.txt"

# R 生物資訊套件
R_PACKAGES="Biostrings GenomicRanges VariantAnnotation"
R_SCRIPT_FILE="bioinformatics-packages.R"

# Conda 環境
CONDA_ENV_FILE="environment-bioinformatics.yml"
```

#### 1.2 更新 Conda 環境檔案
```yaml
# configs/conda/environment-bioinformatics.yml
name: bioinformatics
channels:
  - conda-forge
  - bioconda
  - r
dependencies:
  - python=3.11
  - r-base=4.3
  - fastqc
  - samtools
  - bcftools
  - bedtools
  - biopython
  - pysam
  - pip
  - pip:
    - pyvcf
    - pandas
    - matplotlib
    - seaborn
```

### 2. 專用工具容器 (Phase 2)

#### 2.1 擴展 docker-compose.yml
```yaml
# 在現有 docker-compose.yml 中新增
services:
  # 現有服務...
  devcontainer:
    # 現有配置...
    
  # 新增專用工具服務
  gatk-tools:
    image: broadinstitute/gatk:latest
    volumes:
      - ./data:/workspace/data
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["bioinformatics", "gatk"]
    environment:
      - JAVA_OPTS=-Xmx8g
    
  colabfold:
    image: colabfold/colabfold:latest
    volumes:
      - ./data:/workspace/data
      - ./models:/workspace/models
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["structure", "colabfold"]
    runtime: nvidia  # 需要 GPU 支援
    environment:
      - CUDA_VISIBLE_DEVICES=0
    
  autodock:
    image: ccsb/autodock-vina:latest
    volumes:
      - ./data:/workspace/data
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["molecular", "autodock"]
    
  # 工具管理服務
  tool-manager:
    image: alpine:latest
    volumes:
      - ./:/workspace
      - /var/run/docker.sock:/var/run/docker.sock
    working_dir: /workspace
    profiles: ["tools"]
    command: tail -f /dev/null
```

#### 2.2 工具管理腳本
```bash
#!/bin/bash
# scripts/manage/tool-manager.sh

# GATK 工具包裝器
run_gatk() {
    echo "🔬 執行 GATK: $@"
    docker-compose --profile gatk run --rm gatk-tools gatk "$@"
}

# ColabFold 工具包裝器
run_colabfold() {
    echo "🧬 執行 ColabFold: $@"
    
    # 檢查 GPU 可用性
    if ! nvidia-smi &>/dev/null; then
        echo "⚠️  警告: 未檢測到 GPU，將使用 CPU 模式"
        docker-compose run --rm colabfold colabfold_batch --num-models 1 "$@"
    else
        docker-compose --profile colabfold run --rm colabfold colabfold_batch "$@"
    fi
}

# AutoDock Vina 工具包裝器
run_autodock() {
    echo "💊 執行 AutoDock Vina: $@"
    docker-compose --profile autodock run --rm autodock vina "$@"
}

# 工具狀態檢查
check_tool_availability() {
    local tool="$1"
    
    case $tool in
        "gatk")
            docker-compose --profile gatk run --rm gatk-tools gatk --version
            ;;
        "colabfold")
            docker-compose --profile colabfold run --rm colabfold colabfold_batch --help | head -5
            ;;
        "autodock")
            docker-compose --profile autodock run --rm autodock vina --version
            ;;
        *)
            echo "❌ 未知工具: $tool"
            echo "可用工具: gatk, colabfold, autodock"
            ;;
    esac
}

# 批次工具管理
manage_tools() {
    local action="$1"
    local tools="${2:-all}"
    
    case $action in
        "start")
            echo "🚀 啟動工具容器..."
            docker-compose --profile bioinformatics --profile structure --profile molecular up -d
            ;;
        "stop")
            echo "🛑 停止工具容器..."
            docker-compose --profile bioinformatics --profile structure --profile molecular down
            ;;
        "status")
            echo "📊 工具容器狀態:"
            docker-compose ps
            ;;
        *)
            echo "用法: manage_tools {start|stop|status}"
            ;;
    esac
}

# 主要命令分發
case "$1" in
    "gatk")
        shift
        run_gatk "$@"
        ;;
    "colabfold")
        shift
        run_colabfold "$@"
        ;;
    "autodock")
        shift
        run_autodock "$@"
        ;;
    "check")
        check_tool_availability "$2"
        ;;
    "manage")
        shift
        manage_tools "$@"
        ;;
    *)
        echo "🧬 生物資訊工具管理器"
        echo ""
        echo "輕量級工具 (dev container 內):"
        echo "  fastqc <files>     - 品質控制分析"
        echo "  samtools <cmd>     - SAM/BAM 操作"
        echo "  bcftools <cmd>     - VCF 操作"
        echo "  bedtools <cmd>     - 基因組區間操作"
        echo ""
        echo "重量級工具 (專用容器):"
        echo "  $0 gatk <args>     - GATK 變異呼叫"
        echo "  $0 colabfold <args> - 蛋白質結構預測"
        echo "  $0 autodock <args> - 分子對接"
        echo ""
        echo "管理命令:"
        echo "  $0 check <tool>    - 檢查工具可用性"
        echo "  $0 manage <action> - 管理工具容器"
        ;;
esac
```

### 3. 工作流程整合

#### 3.1 NGS 資料分析工作流程
```bash
#!/bin/bash
# workflows/ngs-analysis.sh

# NGS 資料分析標準流程
ngs_workflow() {
    local input_dir="$1"
    local output_dir="$2"
    local reference="$3"
    
    echo "🧬 開始 NGS 資料分析流程"
    echo "輸入: $input_dir"
    echo "輸出: $output_dir"
    echo "參考基因組: $reference"
    
    # 1. 品質控制 (dev container 內)
    echo "📊 步驟 1: 品質控制"
    mkdir -p "$output_dir/qc"
    fastqc "$input_dir"/*.fastq.gz -o "$output_dir/qc"
    
    # 2. 序列比對 (dev container 內)
    echo "🎯 步驟 2: 序列比對"
    for fastq in "$input_dir"/*.fastq.gz; do
        base=$(basename "$fastq" .fastq.gz)
        bwa mem "$reference" "$fastq" | samtools sort -o "$output_dir/${base}.bam"
        samtools index "$output_dir/${base}.bam"
    done
    
    # 3. 變異呼叫 (GATK 專用容器)
    echo "🔬 步驟 3: 變異呼叫"
    ./scripts/manage/tool-manager.sh gatk HaplotypeCaller \
        -I "$output_dir/sample.bam" \
        -O "$output_dir/variants.vcf" \
        -R "$reference"
    
    # 4. 變異註釋 (dev container 內)
    echo "📝 步驟 4: 變異註釋"
    bcftools annotate "$output_dir/variants.vcf" > "$output_dir/annotated_variants.vcf"
    
    echo "✅ NGS 分析完成!"
}

# 結構預測工作流程
structure_workflow() {
    local protein_fasta="$1"
    local output_dir="$2"
    
    echo "🧬 開始蛋白質結構預測"
    
    # 使用 ColabFold 進行結構預測
    ./scripts/manage/tool-manager.sh colabfold \
        "$protein_fasta" \
        "$output_dir" \
        --num-models 5 \
        --num-recycles 3
    
    echo "✅ 結構預測完成!"
}

# 分子對接工作流程
docking_workflow() {
    local receptor="$1"
    local ligand="$2"
    local output_dir="$3"
    
    echo "💊 開始分子對接分析"
    
    # 使用 AutoDock Vina 進行分子對接
    ./scripts/manage/tool-manager.sh autodock \
        --receptor "$receptor" \
        --ligand "$ligand" \
        --out "$output_dir/docking_result.pdbqt" \
        --log "$output_dir/docking.log"
    
    echo "✅ 分子對接完成!"
}
```

#### 3.2 整合的環境狀態檢查
```bash
#!/bin/bash
# 更新 scripts/manage/devcontainer-status.sh

check_bioinformatics_tools() {
    echo "🧬 生物資訊工具狀態"
    echo "===================="
    
    # 輕量級工具 (dev container 內)
    echo "📊 輕量級工具:"
    echo "  FastQC: $(fastqc --version 2>/dev/null || echo '❌ 未安裝')"
    echo "  SAMtools: $(samtools --version 2>/dev/null | head -1 || echo '❌ 未安裝')"
    echo "  BCFtools: $(bcftools --version 2>/dev/null | head -1 || echo '❌ 未安裝')"
    echo "  BEDtools: $(bedtools --version 2>/dev/null || echo '❌ 未安裝')"
    
    # 重量級工具 (專用容器)
    echo ""
    echo "🔬 重量級工具容器:"
    echo "  GATK: $(docker-compose ps gatk-tools 2>/dev/null | grep -q 'Up' && echo '✅ 運行中' || echo '❌ 未運行')"
    echo "  ColabFold: $(docker-compose ps colabfold 2>/dev/null | grep -q 'Up' && echo '✅ 運行中' || echo '❌ 未運行')"
    echo "  AutoDock: $(docker-compose ps autodock 2>/dev/null | grep -q 'Up' && echo '✅ 運行中' || echo '❌ 未運行')"
    
    # GPU 狀態
    echo ""
    echo "🎮 GPU 狀態:"
    if nvidia-smi &>/dev/null; then
        echo "  NVIDIA GPU: ✅ 可用"
        nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader,nounits | \
            awk -F', ' '{printf "  %s: %dMB/%dMB\n", $1, $3, $2}'
    else
        echo "  NVIDIA GPU: ❌ 不可用 (ColabFold 將使用 CPU)"
    fi
}
```

### 4. 更新的 Profile 系統

#### 4.1 新增生物資訊專用 Profile
```bash
# profiles/bioinformatics-full.env
PROFILE_NAME="Bioinformatics Full"
PROFILE_DESCRIPTION="完整生物資訊分析環境，包含所有工具"

# 輕量級工具
ENABLE_FASTQC=true
ENABLE_SAMTOOLS=true
ENABLE_BCFTOOLS=true
ENABLE_BEDTOOLS=true

# 重量級工具
ENABLE_GATK=true
ENABLE_COLABFOLD=true
ENABLE_AUTODOCK=true

# Docker Compose profiles
COMPOSE_PROFILES="bioinformatics,structure,molecular"

# 資源配置
GATK_MEMORY=8g
COLABFOLD_GPU=true
AUTODOCK_CPU=4
```

#### 4.2 工具選擇性安裝
```bash
# 更新 setup-environment.sh
echo "🧬 生物資訊工具選擇："
echo "請選擇需要的工具組合："
echo ""
echo "1) basic-bio      - 基本工具 (FastQC, SAMtools)"
echo "2) ngs-analysis   - NGS 分析 (basic + GATK)"
echo "3) structure      - 結構預測 (basic + ColabFold)"
echo "4) molecular      - 分子對接 (basic + AutoDock)"
echo "5) full-bio       - 完整生物資訊環境"
echo ""
```

## 📊 資源使用優化

### 1. 按需啟動策略
```bash
# 智能容器管理
smart_container_management() {
    local workflow="$1"
    
    case $workflow in
        "ngs")
            docker-compose --profile bioinformatics up -d
            ;;
        "structure")
            docker-compose --profile structure up -d
            ;;
        "molecular")
            docker-compose --profile molecular up -d
            ;;
        "full")
            docker-compose --profile bioinformatics --profile structure --profile molecular up -d
            ;;
    esac
}
```

### 2. 資源監控
```bash
# 資源使用監控
monitor_resources() {
    echo "📊 容器資源使用情況:"
    docker stats --no-stream --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"
    
    echo ""
    echo "💾 磁碟使用情況:"
    df -h | grep -E "(workspace|models|data)"
}
```

## 🎯 實施優先順序

### Phase 1: 基礎工具整合 (立即可行)
1. ✅ 將 FastQC, SAMtools 整合到現有 bioinformatics profile
2. ✅ 更新 conda 環境檔案
3. ✅ 測試基本工具功能

### Phase 2: 專用容器建置 (中期)
1. 🔄 設定 GATK 專用容器
2. 🔄 建立工具管理腳本
3. 🔄 測試容器間協作

### Phase 3: 進階工具整合 (長期)
1. 🔄 ColabFold GPU 支援設定
2. 🔄 AutoDock 分子對接環境
3. 🔄 完整工作流程整合

### Phase 4: 優化與自動化 (進階)
1. 🔄 資源使用優化
2. 🔄 自動化工作流程
3. 🔄 監控和日誌系統

## 🔧 使用範例

### 日常 NGS 分析流程
```bash
# 1. 啟動 dev container (包含輕量級工具)
code .

# 2. 品質控制 (dev container 內)
fastqc data/*.fastq.gz -o results/qc/

# 3. 基本處理 (dev container 內)
samtools view -h sample.bam | head

# 4. 變異呼叫 (專用容器)
./scripts/manage/tool-manager.sh gatk HaplotypeCaller \
    -I sample.bam -O variants.vcf -R reference.fa

# 5. 後處理 (dev container 內)
bcftools stats variants.vcf > variants_stats.txt
```

### 蛋白質結構預測
```bash
# 1. 檢查 GPU 可用性
./scripts/manage/tool-manager.sh check colabfold

# 2. 執行結構預測
./scripts/manage/tool-manager.sh colabfold \
    protein.fasta \
    results/structure/ \
    --num-models 5
```

## 🏆 預期效益

### 1. 彈性與效率
- **按需使用**: 只啟動需要的工具容器
- **快速切換**: 在不同分析工作流程間快速切換
- **資源優化**: 避免不必要的資源消耗

### 2. 維護性
- **模組化設計**: 每個工具獨立管理
- **版本控制**: 各工具版本獨立更新
- **故障隔離**: 單一工具問題不影響整體環境

### 3. 擴展性
- **新工具整合**: 容易新增新的生物資訊工具
- **客製化工作流程**: 支援不同研究需求
- **多平台支援**: 支援有/無 GPU 的環境

---

**生物資訊工具整合完成後，您將擁有一個既靈活又強大的分析環境，能夠處理從基本的 NGS 分析到複雜的蛋白質結構預測等各種生物資訊任務。**