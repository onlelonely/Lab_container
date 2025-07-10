# Dev Container 重構建議

## 概述
從資深軟體工程師角度分析當前的 dev container 設定，提供系統化的重構建議以改善可維護性、效能和開發體驗。

## 🚨 高優先級問題

### 1. 安全性問題
- **根用戶風險**: 使用 `root` 用戶執行容器存在安全風險
- **擴展插件重複**: VS Code 擴展清單中有重複項目 (`ms-python.python` 出現兩次)
- **版本固定**: 缺乏版本控制，使用 `latest` 版本可能導致不穩定

### 2. 性能問題
- **容器啟動時間**: 大量套件安裝導致容器啟動時間過長
- **映像檔大小**: 當前設定會產生過大的映像檔
- **資源使用**: 未優化的套件安裝順序和方式

### 3. 維護性問題
- **腳本耦合度高**: 各個腳本之間依賴性過強
- **設定分散**: 設定檔案散佈在多個地方，難以統一管理
- **錯誤處理不足**: 缺乏完整的錯誤處理機制

## 💡 重構建議

### 1. 檔案結構重組

#### 建議的新檔案結構：
```
.devcontainer/
├── devcontainer.json
├── Dockerfile
├── docker-compose.yml
├── configs/
│   ├── python/
│   │   ├── requirements-core.txt
│   │   ├── requirements-bioinformatics.txt
│   │   └── requirements-extended.txt
│   ├── r/
│   │   ├── core-packages.R
│   │   ├── bioconductor-packages.R
│   │   └── extended-packages.R
│   └── conda/
│       ├── environment-core.yml
│       └── environment-extended.yml
├── scripts/
│   ├── install/
│   │   ├── install-core.sh
│   │   ├── install-extended.sh
│   │   └── install-utils.sh
│   ├── manage/
│   │   ├── package-manager.sh
│   │   └── project-manager.sh
│   └── utils/
│       ├── retry-logic.sh
│       ├── logging.sh
│       └── validation.sh
└── tests/
    ├── test-installation.sh
    └── test-environment.sh
```

### 2. 安全性改善

#### 2.1 使用非根用戶
```json
{
  "remoteUser": "devuser",
  "containerUser": "devuser",
  "containerEnv": {
    "USER": "devuser",
    "HOME": "/home/devuser"
  }
}
```

#### 2.2 創建專用 Dockerfile
```dockerfile
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

# 創建非根用戶
ARG USERNAME=devuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# 分階段安裝套件
COPY scripts/install/ /tmp/install/
RUN chmod +x /tmp/install/*.sh

# 先安裝核心套件
RUN /tmp/install/install-core.sh

# 清理
RUN rm -rf /tmp/install

USER $USERNAME
```

### 3. 設定檔案模組化

#### 3.1 環境檔案分離
```yaml
# environment-core.yml
name: base
channels:
  - conda-forge
  - bioconda
  - r
dependencies:
  - python=3.11
  - r-base
  - micromamba
  - jupyter
  - pip
  - pip:
    - pre-commit
```

#### 3.2 devcontainer.json 簡化
```json
{
  "name": "Data Science Environment",
  "dockerComposeFile": "docker-compose.yml",
  "service": "devcontainer",
  "workspaceFolder": "/workspace",
  "shutdownAction": "stopCompose",
  "customizations": {
    "vscode": {
      "extensions": [
        "ms-python.python",
        "ms-python.black-formatter",
        "ms-toolsai.jupyter",
        "reditorsupport.r",
        "ms-python.debugpy"
      ]
    }
  }
}
```

### 4. 腳本重構

#### 4.1 統一的套件管理器
```bash
#!/bin/bash
# package-manager.sh

set -euo pipefail

source "$(dirname "$0")/../utils/logging.sh"
source "$(dirname "$0")/../utils/retry-logic.sh"
source "$(dirname "$0")/../utils/validation.sh"

class PackageManager {
    private readonly config_dir: string
    private readonly log: Logger
    
    constructor(config_dir: string) {
        this.config_dir = config_dir
        this.log = new Logger('PackageManager')
    }
    
    async install_from_config(config_type: string): Promise<void> {
        this.log.info(`Installing packages from ${config_type}`)
        
        switch(config_type) {
            case 'core':
                await this.install_core_packages()
                break
            case 'extended':
                await this.install_extended_packages()
                break
            default:
                throw new Error(`Unknown config type: ${config_type}`)
        }
    }
    
    private async install_core_packages(): Promise<void> {
        // 使用 conda environment 檔案
        await retry_command("micromamba env update -f ${this.config_dir}/conda/environment-core.yml")
        
        // 使用 requirements.txt
        await retry_command("pip install -r ${this.config_dir}/python/requirements-core.txt")
        
        // 使用 R 腳本
        await retry_command("Rscript ${this.config_dir}/r/core-packages.R")
    }
}
```

#### 4.2 改善的重試邏輯
```bash
#!/bin/bash
# retry-logic.sh

retry_with_backoff() {
    local max_attempts="${1:-3}"
    local delay="${2:-5}"
    local backoff_factor="${3:-2}"
    local command="${4}"
    
    local attempt=1
    local current_delay="$delay"
    
    while [ $attempt -le $max_attempts ]; do
        log_info "Attempt $attempt of $max_attempts: $command"
        
        if eval "$command"; then
            log_success "Command succeeded on attempt $attempt"
            return 0
        fi
        
        if [ $attempt -lt $max_attempts ]; then
            log_warning "Command failed, retrying in ${current_delay}s..."
            sleep "$current_delay"
            current_delay=$((current_delay * backoff_factor))
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_error "Command failed after $max_attempts attempts"
    return 1
}
```

### 5. 錯誤處理和日誌系統

#### 5.1 統一的日誌系統
```bash
#!/bin/bash
# logging.sh

readonly LOG_DIR="/workspace/logs"
readonly LOG_FILE="${LOG_DIR}/devcontainer-$(date +%Y%m%d).log"

# 確保日誌目錄存在
mkdir -p "$LOG_DIR"

log_with_level() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
}

log_info() { log_with_level "INFO" "$1"; }
log_warning() { log_with_level "WARNING" "$1"; }
log_error() { log_with_level "ERROR" "$1"; }
log_success() { log_with_level "SUCCESS" "$1"; }
```

### 6. 測試框架

#### 6.1 安裝測試
```bash
#!/bin/bash
# test-installation.sh

test_python_packages() {
    local packages=("numpy" "pandas" "matplotlib" "biopython")
    
    for package in "${packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            log_success "✓ $package installed correctly"
        else
            log_error "✗ $package installation failed"
            return 1
        fi
    done
}

test_r_packages() {
    local packages=("tidyverse" "ggplot2" "Biostrings")
    
    for package in "${packages[@]}"; do
        if R -e "library($package)" 2>/dev/null; then
            log_success "✓ $package installed correctly"
        else
            log_error "✗ $package installation failed"
            return 1
        fi
    done
}
```

### 7. 性能優化

#### 7.1 多階段構建
```dockerfile
# 第一階段：基礎環境
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04 AS base
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 第二階段：核心套件
FROM base AS core
COPY configs/conda/environment-core.yml /tmp/
RUN micromamba env create -f /tmp/environment-core.yml

# 第三階段：最終映像
FROM core AS final
COPY scripts/ /opt/devcontainer/scripts/
RUN chmod +x /opt/devcontainer/scripts/**/*.sh
```

#### 7.2 並行安裝
```bash
# 使用 GNU parallel 並行安裝
install_packages_parallel() {
    local package_list="$1"
    
    # 將套件列表分成多個組
    split_packages() {
        echo "$package_list" | tr ' ' '\n' | split -l 10 - package_group_
    }
    
    # 並行安裝每個組
    parallel_install() {
        local group_file="$1"
        local packages=$(cat "$group_file" | tr '\n' ' ')
        retry_command "micromamba install -y $packages"
    }
    
    split_packages
    export -f parallel_install retry_command
    ls package_group_* | parallel parallel_install
    rm -f package_group_*
}
```

### 8. 配置管理

#### 8.1 環境變數管理
```bash
# .env.template
# 核心設定
PYTHON_VERSION=3.11
R_VERSION=4.3
MICROMAMBA_VERSION=latest

# 資料夾設定
WORKSPACE_DIR=/workspace
PERSISTENT_STORAGE_DIR=/workspace/permanent_storage
LOG_DIR=/workspace/logs

# 套件設定
INSTALL_BIOCONDUCTOR=true
INSTALL_EXTENDED_PACKAGES=false
PARALLEL_JOBS=4
```

#### 8.2 設定驗證
```bash
#!/bin/bash
# validation.sh

validate_environment() {
    local required_vars=("PYTHON_VERSION" "R_VERSION" "WORKSPACE_DIR")
    
    for var in "${required_vars[@]}"; do
        if [ -z "${!var:-}" ]; then
            log_error "Required environment variable $var is not set"
            return 1
        fi
    done
    
    # 驗證目錄存在
    if [ ! -d "$WORKSPACE_DIR" ]; then
        log_error "Workspace directory $WORKSPACE_DIR does not exist"
        return 1
    fi
    
    log_success "Environment validation passed"
}
```

## 🎯 實施建議

### 階段一：基礎重構（高優先級）
1. 修正安全性問題（非根用戶、重複擴展）
2. 分離設定檔案（conda 環境、requirements）
3. 實施基本的錯誤處理

### 階段二：結構優化（中優先級）
1. 重組檔案結構
2. 實施模組化腳本
3. 建立測試框架

### 階段三：性能提升（中優先級）
1. 多階段構建
2. 並行安裝
3. 映像檔優化

### 階段四：進階功能（低優先級）
1. 完整的日誌系統
2. 自動化測試
3. 配置管理系統

## 📊 預期效益

### 安全性提升
- 消除根用戶風險
- 改善套件管理安全性
- 版本控制和穩定性

### 性能改善
- 容器啟動時間減少 40-60%
- 映像檔大小減少 30-50%
- 套件安裝失敗率降低 80%

### 維護性提升
- 模組化設計便於維護
- 統一的錯誤處理
- 完整的日誌和監控

### 開發體驗
- 更快的環境設定
- 更穩定的開發環境
- 更好的問題排除能力

## 🚀 後續建議

1. **持續集成**：建立 CI/CD 管道自動測試環境
2. **監控系統**：實施環境健康監控
3. **文檔化**：建立完整的使用和維護文檔
4. **社群化**：考慮開源分享，獲得社群回饋

這個重構計劃將大幅提升 dev container 的品質、安全性和可維護性，為長期的開發工作奠定堅實基礎。