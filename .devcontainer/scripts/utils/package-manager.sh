#!/bin/bash
# package-manager.sh - 動態套件管理工具

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

# 確保必要目錄存在
ensure_directories() {
    mkdir -p .devcontainer/configs/python
    mkdir -p .devcontainer/configs/r
    mkdir -p .devcontainer/configs/conda
}

# 新增 Python 套件
add_python_package() {
    local package_name="$1"
    local runtime_file=".devcontainer/configs/python/requirements-runtime.txt"
    
    info_msg "安裝 Python 套件: $package_name"
    
    # 安裝套件
    if pip install "$package_name"; then
        # 記錄到運行時檔案
        echo "$package_name" >> "$runtime_file"
        success_msg "Python 套件 '$package_name' 安裝成功"
    else
        error_exit "Python 套件 '$package_name' 安裝失敗"
    fi
}

# 新增 R 套件
add_r_package() {
    local package_name="$1"
    local runtime_file=".devcontainer/configs/r/runtime-packages.R"
    
    info_msg "安裝 R 套件: $package_name"
    
    # 安裝套件
    if R -e "install.packages('$package_name', repos='https://cran.rstudio.com/')"; then
        # 記錄到運行時檔案
        echo "install.packages('$package_name')" >> "$runtime_file"
        success_msg "R 套件 '$package_name' 安裝成功"
    else
        error_exit "R 套件 '$package_name' 安裝失敗"
    fi
}

# 新增 Conda 套件
add_conda_package() {
    local package_name="$1"
    local runtime_file=".devcontainer/configs/conda/runtime-packages.yml"
    
    info_msg "安裝 Conda 套件: $package_name"
    
    # 安裝套件
    if conda install -y "$package_name"; then
        # 記錄到運行時檔案
        if [ ! -f "$runtime_file" ]; then
            cat > "$runtime_file" << EOF
# Runtime installed conda packages
dependencies:
EOF
        fi
        echo "  - $package_name" >> "$runtime_file"
        success_msg "Conda 套件 '$package_name' 安裝成功"
    else
        error_exit "Conda 套件 '$package_name' 安裝失敗"
    fi
}

# 移除 Python 套件
remove_python_package() {
    local package_name="$1"
    local runtime_file=".devcontainer/configs/python/requirements-runtime.txt"
    
    info_msg "移除 Python 套件: $package_name"
    
    # 移除套件
    if pip uninstall -y "$package_name"; then
        # 從運行時檔案移除
        if [ -f "$runtime_file" ]; then
            grep -v "^$package_name$" "$runtime_file" > "$runtime_file.tmp" || true
            mv "$runtime_file.tmp" "$runtime_file"
        fi
        success_msg "Python 套件 '$package_name' 移除成功"
    else
        error_exit "Python 套件 '$package_name' 移除失敗"
    fi
}

# 列出已安裝套件
list_packages() {
    echo ""
    echo -e "${BLUE}📦 已安裝套件清單${NC}"
    echo "===================="
    
    echo ""
    echo -e "${GREEN}🐍 Python 套件${NC}"
    echo "----------------"
    pip list --format=columns
    
    echo ""
    echo -e "${GREEN}📊 R 套件${NC}"
    echo "----------------"
    R -e "installed.packages()[,c('Package','Version')] |> as.data.frame() |> head(20)" --quiet --no-restore
    
    echo ""
    echo -e "${GREEN}🐍 Conda 套件${NC}"
    echo "----------------"
    conda list | head -20
    
    echo ""
    info_msg "顯示前 20 個套件，完整清單請使用 'pip list', 'conda list' 等命令"
}

# 搜尋套件
search_package() {
    local package_name="$1"
    local package_type="${2:-all}"
    
    echo ""
    echo -e "${BLUE}🔍 搜尋套件: $package_name${NC}"
    echo "=========================="
    
    case $package_type in
        "python"|"all")
            echo ""
            echo -e "${GREEN}🐍 Python 套件搜尋結果${NC}"
            pip search "$package_name" 2>/dev/null || echo "PyPI 搜尋服務暫時不可用，請直接嘗試安裝"
            ;;
        "conda"|"all")
            echo ""
            echo -e "${GREEN}🐍 Conda 套件搜尋結果${NC}"
            conda search "$package_name" 2>/dev/null || echo "找不到 conda 套件"
            ;;
    esac
}

# 備份當前環境
backup_environment() {
    local backup_dir="backups/$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$backup_dir"
    
    info_msg "備份環境到: $backup_dir"
    
    # 備份 Python 套件
    pip freeze > "$backup_dir/requirements.txt"
    
    # 備份 R 套件
    R -e "installed.packages()[,c('Package','Version')] |> write.csv('$backup_dir/r-packages.csv')" --quiet --no-restore
    
    # 備份 Conda 環境
    conda env export > "$backup_dir/environment.yml"
    
    # 備份當前 profile
    if [ -f ".devcontainer/.current-profile" ]; then
        cp ".devcontainer/.current-profile" "$backup_dir/"
    fi
    
    success_msg "環境已備份至: $backup_dir"
}

# 還原環境
restore_environment() {
    local backup_dir="$1"
    
    if [ ! -d "$backup_dir" ]; then
        error_exit "備份目錄不存在: $backup_dir"
    fi
    
    warning_msg "這將會覆蓋當前環境，確定要繼續嗎？"
    read -p "請輸入 'yes' 確認: " confirm
    
    if [ "$confirm" != "yes" ]; then
        info_msg "操作已取消"
        return 0
    fi
    
    info_msg "從備份還原環境: $backup_dir"
    
    # 還原 Python 套件
    if [ -f "$backup_dir/requirements.txt" ]; then
        pip install -r "$backup_dir/requirements.txt"
    fi
    
    # 還原 Conda 環境
    if [ -f "$backup_dir/environment.yml" ]; then
        conda env update --file "$backup_dir/environment.yml"
    fi
    
    # 還原 profile
    if [ -f "$backup_dir/.current-profile" ]; then
        cp "$backup_dir/.current-profile" ".devcontainer/"
    fi
    
    success_msg "環境還原完成"
}

# 顯示使用說明
show_help() {
    echo ""
    echo -e "${BLUE}📦 動態套件管理工具${NC}"
    echo "======================="
    echo ""
    echo "新增套件:"
    echo "  $0 add python <package>  - 新增 Python 套件"
    echo "  $0 add r <package>       - 新增 R 套件"
    echo "  $0 add conda <package>   - 新增 Conda 套件"
    echo ""
    echo "移除套件:"
    echo "  $0 remove python <package> - 移除 Python 套件"
    echo ""
    echo "管理功能:"
    echo "  $0 list                  - 列出已安裝套件"
    echo "  $0 search <package>      - 搜尋套件"
    echo "  $0 backup                - 備份當前環境"
    echo "  $0 restore <backup_dir>  - 還原環境"
    echo ""
    echo "範例:"
    echo "  $0 add python requests"
    echo "  $0 add r forecast"
    echo "  $0 add conda nodejs"
    echo "  $0 search tensorflow"
    echo ""
}

# 主程式
main() {
    # 確保目錄存在
    ensure_directories
    
    # 檢查參數
    if [ $# -eq 0 ]; then
        show_help
        return 0
    fi
    
    local action="$1"
    
    case $action in
        "add")
            if [ $# -lt 3 ]; then
                error_exit "用法: $0 add <type> <package>"
            fi
            local package_type="$2"
            local package_name="$3"
            
            case $package_type in
                "python") add_python_package "$package_name" ;;
                "r") add_r_package "$package_name" ;;
                "conda") add_conda_package "$package_name" ;;
                *) error_exit "不支援的套件類型: $package_type" ;;
            esac
            ;;
        "remove")
            if [ $# -lt 3 ]; then
                error_exit "用法: $0 remove <type> <package>"
            fi
            local package_type="$2"
            local package_name="$3"
            
            case $package_type in
                "python") remove_python_package "$package_name" ;;
                *) error_exit "目前只支援移除 Python 套件" ;;
            esac
            ;;
        "list")
            list_packages
            ;;
        "search")
            if [ $# -lt 2 ]; then
                error_exit "用法: $0 search <package> [type]"
            fi
            local package_name="$2"
            local package_type="${3:-all}"
            search_package "$package_name" "$package_type"
            ;;
        "backup")
            backup_environment
            ;;
        "restore")
            if [ $# -lt 2 ]; then
                error_exit "用法: $0 restore <backup_dir>"
            fi
            local backup_dir="$2"
            restore_environment "$backup_dir"
            ;;
        "help"|"-h"|"--help")
            show_help
            ;;
        *)
            error_exit "未知動作: $action"
            ;;
    esac
}

# 執行主程式
main "$@"