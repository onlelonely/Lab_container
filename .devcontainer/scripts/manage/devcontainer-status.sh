#!/bin/bash
# devcontainer-status.sh - Dev Container 環境狀態儀表板

set -euo pipefail

# 顏色定義
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# 取得當前 profile 資訊
get_current_profile() {
    if [ -f ".devcontainer/.current-profile" ]; then
        local current=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current}.env"
        
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            echo "$PROFILE_NAME ($current)"
        else
            echo "未知 ($current)"
        fi
    else
        echo "未設定"
    fi
}

# 檢查 Python 版本和套件
check_python() {
    echo -e "${GREEN}🐍 Python 環境${NC}"
    echo "----------------"
    
    if command -v python >/dev/null 2>&1; then
        echo "版本: $(python --version)"
        echo "位置: $(which python)"
        local package_count=$(pip list 2>/dev/null | wc -l)
        echo "已安裝套件: $package_count 個"
        
        # 檢查常用套件
        echo "核心套件狀態:"
        for pkg in pandas numpy matplotlib jupyter; do
            if pip show "$pkg" >/dev/null 2>&1; then
                local version=$(pip show "$pkg" | grep Version | cut -d' ' -f2)
                echo "  ✅ $pkg ($version)"
            else
                echo "  ❌ $pkg (未安裝)"
            fi
        done
    else
        echo "❌ Python 未安裝"
    fi
}

# 檢查 R 版本和套件
check_r() {
    echo ""
    echo -e "${GREEN}📊 R 環境${NC}"
    echo "----------------"
    
    if command -v R >/dev/null 2>&1; then
        echo "版本: $(R --version | head -1 | cut -d' ' -f3-4)"
        echo "位置: $(which R)"
        
        # 檢查 R 套件數量
        local r_packages=$(R -e "length(installed.packages()[,1])" --quiet --slave 2>/dev/null | tail -1)
        echo "已安裝套件: $r_packages 個"
        
        # 檢查常用 R 套件
        echo "核心套件狀態:"
        for pkg in "tidyverse" "ggplot2" "dplyr"; do
            if R -e "library($pkg)" --quiet --slave >/dev/null 2>&1; then
                echo "  ✅ $pkg"
            else
                echo "  ❌ $pkg (未安裝)"
            fi
        done
    else
        echo "❌ R 未安裝"
    fi
}

# 檢查 Jupyter 環境
check_jupyter() {
    echo ""
    echo -e "${GREEN}📚 Jupyter 環境${NC}"
    echo "----------------"
    
    if command -v jupyter >/dev/null 2>&1; then
        echo "JupyterLab: $(jupyter --version 2>/dev/null | grep jupyterlab || echo '未安裝')"
        echo "Notebook: $(jupyter --version 2>/dev/null | grep notebook || echo '未安裝')"
        
        # 檢查 kernels
        echo "可用 Kernels:"
        jupyter kernelspec list 2>/dev/null | grep -v "Available kernels:" | while read line; do
            if [ -n "$line" ]; then
                echo "  📋 $line"
            fi
        done
    else
        echo "❌ Jupyter 未安裝"
    fi
}

# 檢查生物資訊工具
check_bioinformatics_tools() {
    echo ""
    echo -e "${GREEN}🧬 生物資訊工具${NC}"
    echo "----------------"
    
    # 檢查輕量級工具
    echo "輕量級工具:"
    for tool in fastqc samtools bcftools bedtools; do
        if command -v "$tool" >/dev/null 2>&1; then
            local version=$($tool --version 2>/dev/null | head -1 || echo "unknown")
            echo "  ✅ $tool ($version)"
        else
            echo "  ❌ $tool (未安裝)"
        fi
    done
    
    # 檢查 Python 生物資訊套件
    echo "Python 生物資訊套件:"
    for pkg in biopython pysam; do
        if pip show "$pkg" >/dev/null 2>&1; then
            local version=$(pip show "$pkg" | grep Version | cut -d' ' -f2)
            echo "  ✅ $pkg ($version)"
        else
            echo "  ❌ $pkg (未安裝)"
        fi
    done
    
    # 檢查 R 生物資訊套件
    echo "R 生物資訊套件:"
    for pkg in "Biostrings" "GenomicRanges"; do
        if R -e "library($pkg)" --quiet --slave >/dev/null 2>&1; then
            echo "  ✅ $pkg"
        else
            echo "  ❌ $pkg (未安裝)"
        fi
    done
}

# 檢查系統資源
check_system_resources() {
    echo ""
    echo -e "${GREEN}💻 系統資源${NC}"
    echo "----------------"
    
    # CPU 資訊
    echo "CPU: $(nproc) 核心"
    
    # 記憶體資訊
    if command -v free >/dev/null 2>&1; then
        local mem_info=$(free -h | grep Mem)
        local total=$(echo $mem_info | awk '{print $2}')
        local used=$(echo $mem_info | awk '{print $3}')
        local available=$(echo $mem_info | awk '{print $7}')
        echo "記憶體: $used / $total (可用: $available)"
    fi
    
    # 磁碟空間
    echo "磁碟使用:"
    if [ -d "/workspace" ]; then
        df -h /workspace 2>/dev/null | tail -1 | awk '{printf "  /workspace: %s / %s (%s 已使用)\n", $3, $2, $5}'
    else
        df -h . 2>/dev/null | tail -1 | awk '{printf "  當前目錄: %s / %s (%s 已使用)\n", $3, $2, $5}'
    fi
}

# 檢查網路連線
check_network() {
    echo ""
    echo -e "${GREEN}🌐 網路連線${NC}"
    echo "----------------"
    
    # 檢查基本網路連線
    if ping -c 1 8.8.8.8 >/dev/null 2>&1; then
        echo "網路連線: ✅ 正常"
    else
        echo "網路連線: ❌ 異常"
    fi
    
    # 檢查常用套件庫連線
    echo "套件庫連線:"
    
    # PyPI
    if curl -s --max-time 3 https://pypi.org >/dev/null 2>&1; then
        echo "  ✅ PyPI (Python)"
    else
        echo "  ❌ PyPI (Python)"
    fi
    
    # CRAN
    if curl -s --max-time 3 https://cran.r-project.org >/dev/null 2>&1; then
        echo "  ✅ CRAN (R)"
    else
        echo "  ❌ CRAN (R)"
    fi
    
    # Conda
    if curl -s --max-time 3 https://conda.anaconda.org >/dev/null 2>&1; then
        echo "  ✅ Conda"
    else
        echo "  ❌ Conda"
    fi
}

# 檢查專用工具容器狀態
check_tool_containers() {
    echo ""
    echo -e "${GREEN}🔬 專用工具容器${NC}"
    echo "----------------"
    
    # 檢查是否有 docker-compose.tools.yml
    if [ -f ".devcontainer/docker-compose.tools.yml" ]; then
        echo "工具容器配置: ✅ 可用"
        
        # 檢查各工具容器狀態
        local tools=("gatk-tools" "colabfold" "autodock")
        for tool in "${tools[@]}"; do
            if docker ps --format "table {{.Names}}" | grep -q "devcontainer-${tool%-*}"; then
                echo "  ✅ $tool (運行中)"
            else
                echo "  ❌ $tool (未運行)"
            fi
        done
        
        # 檢查工具管理器
        if [ -f ".devcontainer/scripts/manage/tool-manager.sh" ]; then
            echo "工具管理器: ✅ 可用"
        else
            echo "工具管理器: ❌ 不可用"
        fi
    else
        echo "工具容器配置: ❌ 未設定"
    fi
}

# 檢查容器狀態
check_container_status() {
    echo ""
    echo -e "${GREEN}🐳 容器資訊${NC}"
    echo "----------------"
    
    # 檢查是否在容器內
    if [ -f /.dockerenv ] || [ -n "${DEVCONTAINER:-}" ]; then
        echo "環境: 🐳 Dev Container"
        
        # 容器 ID
        if [ -f /.dockerenv ]; then
            local container_id=$(basename $(cat /proc/1/cpuset) | head -c 12)
            echo "容器 ID: $container_id"
        fi
        
        # 檢查掛載點
        echo "掛載點:"
        mount | grep -E "(workspace|home)" | while read line; do
            echo "  📁 $line"
        done | head -3
        
    else
        echo "環境: 💻 本機"
    fi
}

# 顯示快速命令
show_quick_commands() {
    echo ""
    echo -e "${CYAN}⚡ 快速命令${NC}"
    echo "----------------"
    echo "環境管理:"
    echo "  bash .devcontainer/scripts/manage/setup-environment.sh"
    echo ""
    echo "套件管理:"
    echo "  bash .devcontainer/scripts/utils/package-manager.sh"
    echo ""
    echo "工具管理:"
    echo "  bash .devcontainer/scripts/manage/tool-manager.sh"
    echo ""
    echo "檢查狀態:"
    echo "  bash .devcontainer/scripts/manage/devcontainer-status.sh"
    echo ""
}

# 顯示運行時資訊
show_runtime_info() {
    echo ""
    echo -e "${YELLOW}⏰ 運行時資訊${NC}"
    echo "----------------"
    echo "當前時間: $(date)"
    echo "系統啟動: $(uptime -s 2>/dev/null || echo '無法取得')"
    echo "負載平均: $(uptime | awk -F'load average:' '{print $2}' | xargs)"
}

# 主要狀態檢查函數
main_status() {
    # 標題
    echo ""
    echo -e "${BLUE}📊 Dev Container 環境狀態${NC}"
    echo "=============================================="
    
    # 當前 profile
    echo ""
    echo -e "${MAGENTA}🔧 當前環境設定${NC}"
    echo "----------------"
    echo "Profile: $(get_current_profile)"
    
    # 各項檢查
    check_python
    check_r
    check_jupyter
    
    # 如果當前 profile 啟用了生物資訊工具，則檢查
    if [ -f ".devcontainer/.current-profile" ]; then
        local current=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current}.env"
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
                check_bioinformatics_tools
            fi
        fi
    fi
    
    check_system_resources
    check_network
    check_container_status
    check_tool_containers
    show_runtime_info
    show_quick_commands
    
    echo ""
    echo -e "${GREEN}✅ 狀態檢查完成${NC}"
}

# 簡化檢查
quick_status() {
    echo -e "${BLUE}📊 快速狀態檢查${NC}"
    echo "==================="
    echo "Profile: $(get_current_profile)"
    echo "Python: $(python --version 2>/dev/null || echo '❌ 未安裝')"
    echo "R: $(R --version 2>/dev/null | head -1 | cut -d' ' -f3-4 || echo '❌ 未安裝')"
    echo "Jupyter: $(jupyter --version 2>/dev/null | head -1 || echo '❌ 未安裝')"
    echo "磁碟: $(df -h . | tail -1 | awk '{print $5 " 已使用"}')"
    echo "記憶體: $(free | grep Mem | awk '{printf "%.1f%% 已使用\n", $3/$2 * 100.0}' 2>/dev/null || echo '無法取得')"
}

# 主程式
main() {
    local mode="${1:-full}"
    
    case $mode in
        "quick"|"-q"|"--quick")
            quick_status
            ;;
        "full"|"")
            main_status
            ;;
        "help"|"-h"|"--help")
            echo "用法: $0 [quick|full]"
            echo "  quick - 快速狀態檢查"
            echo "  full  - 完整狀態檢查 (預設)"
            ;;
        *)
            echo "未知選項: $mode"
            echo "用法: $0 [quick|full]"
            exit 1
            ;;
    esac
}

# 執行主程式
main "$@"