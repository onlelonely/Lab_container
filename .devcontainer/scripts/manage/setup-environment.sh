#!/bin/bash
# setup-environment.sh - 互動式環境選擇器

set -euo pipefail

# 顏色定義
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 錯誤處理
error_exit() {
    echo -e "${RED}❌ 錯誤: $1${NC}" >&2
    exit 1
}

# 成功訊息
success_msg() {
    echo -e "${GREEN}✅ $1${NC}"
}

# 警告訊息
warning_msg() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

# 資訊訊息
info_msg() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

# 載入 profile 函數
load_profile() {
    local profile_name="$1"
    local profile_file=".devcontainer/profiles/${profile_name}.env"
    
    if [ ! -f "$profile_file" ]; then
        error_exit "找不到環境設定檔: $profile_file"
    fi
    
    # 載入環境變數
    set -a
    source "$profile_file"
    set +a
    
    success_msg "載入環境設定: $PROFILE_NAME"
    
    # 儲存當前選擇的 profile
    echo "$profile_name" > .devcontainer/.current-profile
    
    # 生成配置檔案
    generate_configurations
}

# 生成配置檔案
generate_configurations() {
    info_msg "生成配置檔案..."
    
    # 生成 requirements.txt (如果指定了)
    if [ -n "${PYTHON_REQUIREMENTS_FILE:-}" ]; then
        generate_python_requirements
    fi
    
    # 生成 R packages 安裝腳本
    if [ -n "${R_SCRIPT_FILE:-}" ]; then
        generate_r_packages
    fi
    
    # 生成 conda 環境檔案
    if [ -n "${CONDA_ENV_FILE:-}" ]; then
        generate_conda_environment
    fi
    
    # 生成 Docker 配置範本
    generate_docker_templates
    
    success_msg "配置檔案生成完成"
}

# 生成 Python requirements
generate_python_requirements() {
    local req_file=".devcontainer/configs/python/${PYTHON_REQUIREMENTS_FILE}"
    mkdir -p "$(dirname "$req_file")"
    
    echo "# Auto-generated requirements for profile: $PROFILE_NAME" > "$req_file"
    echo "# Generated on: $(date)" >> "$req_file"
    echo "" >> "$req_file"
    
    # 將套件列表轉換為 requirements 格式
    for package in $PYTHON_PACKAGES; do
        echo "$package" >> "$req_file"
    done
    
    info_msg "已生成 Python requirements: $req_file"
}

# 生成 R packages 安裝腳本
generate_r_packages() {
    local r_file=".devcontainer/configs/r/${R_SCRIPT_FILE}"
    mkdir -p "$(dirname "$r_file")"
    
    cat > "$r_file" << EOF
# Auto-generated R packages installation script for profile: $PROFILE_NAME
# Generated on: $(date)

# 設定 CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 安裝套件
packages <- c("$( echo "$R_PACKAGES" | sed 's/ /", "/g' )")

for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

cat("✅ R packages installation completed for profile: $PROFILE_NAME\\n")
EOF
    
    info_msg "已生成 R packages 腳本: $r_file"
}

# 生成 conda 環境檔案
generate_conda_environment() {
    local conda_file=".devcontainer/configs/conda/${CONDA_ENV_FILE}"
    mkdir -p "$(dirname "$conda_file")"
    
    cat > "$conda_file" << EOF
# Auto-generated conda environment for profile: $PROFILE_NAME
# Generated on: $(date)
name: ${PROFILE_NAME,,}  # 轉換為小寫
channels:
  - conda-forge
  - bioconda
  - r
dependencies:
  - python=3.11
  - r-base=4.3
EOF

    # 如果啟用生物資訊工具，添加相關套件
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        cat >> "$conda_file" << EOF
  # 生物資訊工具
EOF
        for tool in ${BIOINFORMATICS_TOOLS:-}; do
            echo "  - $tool" >> "$conda_file"
        done
    fi
    
    cat >> "$conda_file" << EOF
  - pip
  - pip:
EOF
    
    # 添加 Python 套件
    for package in $PYTHON_PACKAGES; do
        echo "    - $package" >> "$conda_file"
    done
    
    info_msg "已生成 conda 環境檔案: $conda_file"
}

# 生成 Docker 配置範本
generate_docker_templates() {
    info_msg "生成 Docker 配置範本..."
    
    # 執行範本生成器
    if [ -f ".devcontainer/scripts/utils/template-generator.sh" ]; then
        bash .devcontainer/scripts/utils/template-generator.sh
    else
        warning_msg "找不到範本生成器，跳過 Docker 配置生成"
    fi
}

# 自訂環境設定
setup_custom_environment() {
    echo ""
    echo -e "${BLUE}🔧 自訂環境設定${NC}"
    echo "請回答以下問題來建立您的客製化環境："
    echo ""
    
    read -p "環境名稱: " custom_name
    read -p "環境描述: " custom_desc
    
    echo ""
    echo "請選擇需要的 Python 套件 (用空格分隔):"
    echo "建議: pandas numpy matplotlib seaborn scikit-learn"
    read -p "Python 套件: " custom_python
    
    echo ""
    echo "請選擇需要的 R 套件 (用空格分隔):"
    echo "建議: tidyverse ggplot2 dplyr"
    read -p "R 套件: " custom_r
    
    echo ""
    read -p "是否需要生物資訊工具? (y/N): " need_bio
    
    # 建立自訂 profile
    local custom_file=".devcontainer/profiles/custom.env"
    cat > "$custom_file" << EOF
# custom.env - 自訂環境
PROFILE_NAME="$custom_name"
PROFILE_DESCRIPTION="$custom_desc"

# Python 套件
PYTHON_PACKAGES="$custom_python"
PYTHON_REQUIREMENTS_FILE="requirements-custom.txt"

# R 套件
R_PACKAGES="$custom_r"
R_SCRIPT_FILE="custom-packages.R"

# Conda 環境
CONDA_ENV_FILE="environment-custom.yml"

# 啟用的功能
ENABLE_JUPYTER=true
ENABLE_RSTUDIO=true
ENABLE_VSCODE_R=true
ENABLE_VSCODE_PYTHON=true

# 生物資訊工具
EOF
    
    if [[ "$need_bio" =~ ^[Yy]$ ]]; then
        cat >> "$custom_file" << EOF
ENABLE_BIOINFORMATICS=true
BIOINFORMATICS_TOOLS="fastqc samtools bcftools bedtools"
EOF
    else
        cat >> "$custom_file" << EOF
ENABLE_BIOINFORMATICS=false
EOF
    fi
    
    success_msg "自訂環境設定已建立"
    load_profile "custom"
}

# 顯示當前環境資訊
show_current_environment() {
    if [ -f ".devcontainer/.current-profile" ]; then
        local current=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current}.env"
        
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            echo ""
            echo -e "${GREEN}📊 當前環境資訊${NC}"
            echo "=========================="
            echo "環境名稱: $PROFILE_NAME"
            echo "環境描述: $PROFILE_DESCRIPTION"
            echo "Profile 檔案: $current"
            echo ""
        fi
    else
        warning_msg "尚未設定環境"
    fi
}

# 主選單
main_menu() {
    echo ""
    echo -e "${BLUE}🔧 Dev Container 環境設定${NC}"
    echo "請選擇您需要的分析環境："
    echo ""
    echo "1) minimal        - 基本 Python/R 環境"
    echo "2) datascience    - 資料科學標準環境 (pandas, tidyverse)"
    echo "3) bioinformatics - NGS 生物資訊環境 (Bioconductor, BioPython)"
    echo "4) ml             - 機器學習環境 (scikit-learn, tensorflow)"
    echo "5) statistics     - 統計分析環境 (R focus)"
    echo "6) full           - 完整功能環境"
    echo "7) custom         - 自訂環境"
    echo "8) show           - 顯示當前環境"
    echo "9) exit           - 離開"
    echo ""
    
    read -p "請輸入選項 (1-9): " choice
    
    case $choice in
        1) load_profile "minimal" ;;
        2) load_profile "datascience" ;;
        3) load_profile "bioinformatics" ;;
        4) load_profile "ml" ;;
        5) load_profile "statistics" ;;
        6) load_profile "full" ;;
        7) setup_custom_environment ;;
        8) show_current_environment && main_menu ;;
        9) 
            echo "再見！"
            exit 0
            ;;
        *) 
            warning_msg "無效選項，請重新選擇"
            main_menu
            ;;
    esac
}

# 檢查必要目錄
check_directories() {
    mkdir -p .devcontainer/configs/python
    mkdir -p .devcontainer/configs/r
    mkdir -p .devcontainer/configs/conda
}

# 主程式
main() {
    # 檢查是否在正確的目錄
    if [ ! -d ".devcontainer" ]; then
        error_exit "請在包含 .devcontainer 目錄的專案根目錄執行此腳本"
    fi
    
    # 建立必要目錄
    check_directories
    
    # 如果有命令列參數，直接載入指定的 profile
    if [ $# -gt 0 ]; then
        local profile_name="$1"
        if [ -f ".devcontainer/profiles/${profile_name}.env" ]; then
            load_profile "$profile_name"
            success_msg "環境設定完成！"
            return 0
        else
            error_exit "找不到 profile: $profile_name"
        fi
    fi
    
    # 顯示歡迎訊息
    echo -e "${GREEN}🚀 歡迎使用 Dev Container 環境管理工具${NC}"
    
    # 顯示當前環境（如果有的話）
    show_current_environment
    
    # 顯示主選單
    main_menu
    
    success_msg "環境設定完成！"
    info_msg "您可以現在啟動 VS Code 或重建 dev container 來使用新環境"
}

# 執行主程式
main "$@"