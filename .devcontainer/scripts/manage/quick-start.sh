#!/bin/bash
# quick-start.sh - Dev Container 快速啟動腳本

set -euo pipefail

# 顏色定義
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 顯示歡迎訊息
show_welcome() {
    echo ""
    echo -e "${CYAN}🚀 Dev Container 快速啟動工具${NC}"
    echo "======================================="
    echo ""
    echo -e "${GREEN}歡迎使用彈性分析環境！${NC}"
    echo ""
}

# 檢查現有環境
check_existing_environment() {
    if [ -f ".devcontainer/.current-profile" ]; then
        local current_profile=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current_profile}.env"
        
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            echo -e "${BLUE}📋 發現現有環境配置${NC}"
            echo "環境名稱: $PROFILE_NAME"
            echo "Profile: $current_profile"
            echo "描述: $PROFILE_DESCRIPTION"
            echo ""
            
            read -p "是否要重新設定環境? (y/N): " reset_env
            
            if [[ ! "$reset_env" =~ ^[Yy]$ ]]; then
                echo -e "${GREEN}✅ 使用現有環境設定${NC}"
                return 0
            else
                echo -e "${YELLOW}🔄 重新設定環境...${NC}"
                return 1
            fi
        fi
    fi
    
    echo -e "${YELLOW}🆕 未發現現有環境配置，開始初始設定...${NC}"
    return 1
}

# 執行環境設定
setup_environment() {
    echo -e "${BLUE}🔧 執行環境設定腳本...${NC}"
    bash .devcontainer/scripts/manage/setup-environment.sh
}

# 顯示環境狀態
show_status() {
    echo ""
    echo -e "${CYAN}📊 環境狀態檢查${NC}"
    echo "=================="
    bash .devcontainer/scripts/manage/devcontainer-status.sh quick
}

# 顯示可用工具
show_available_tools() {
    echo ""
    echo -e "${GREEN}🛠️  可用工具和命令${NC}"
    echo "===================="
    echo ""
    echo "📋 環境管理:"
    echo "  setup-env      - 重新設定環境"
    echo "  status         - 檢查環境狀態"
    echo "  quick-start    - 快速啟動"
    echo ""
    echo "📦 套件管理:"
    echo "  add-package    - 新增套件"
    echo "  list-packages  - 列出套件"
    echo "  backup-env     - 備份環境"
    echo ""
    echo "🧬 分析工具 (如果已啟用):"
    
    # 檢查是否啟用生物資訊工具
    if [ -f ".devcontainer/.current-profile" ]; then
        local current=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current}.env"
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
                echo "  fastqc         - 序列品質控制"
                echo "  samtools       - SAM/BAM 操作"
                echo "  bcftools       - VCF 操作"
                echo "  bedtools       - 基因組區間操作"
            fi
        fi
    fi
    
    echo ""
    echo "💡 完整說明請執行: bash .devcontainer/scripts/manage/setup-environment.sh"
}

# 建立快捷命令別名
create_aliases() {
    local alias_file="$HOME/.devcontainer_aliases"
    
    cat > "$alias_file" << 'EOF'
# Dev Container 快捷命令別名

# 環境管理
alias setup-env='bash .devcontainer/scripts/manage/setup-environment.sh'
alias status='bash .devcontainer/scripts/manage/devcontainer-status.sh'
alias quick-start='bash .devcontainer/scripts/manage/quick-start.sh'

# 套件管理
alias add-package='bash .devcontainer/scripts/utils/package-manager.sh add'
alias list-packages='bash .devcontainer/scripts/utils/package-manager.sh list'
alias backup-env='bash .devcontainer/scripts/utils/package-manager.sh backup'

# 快速狀態檢查
alias qs='bash .devcontainer/scripts/manage/devcontainer-status.sh quick'
EOF
    
    # 將別名添加到 bashrc
    if ! grep -q "source.*devcontainer_aliases" "$HOME/.bashrc" 2>/dev/null; then
        echo "" >> "$HOME/.bashrc"
        echo "# Dev Container 別名" >> "$HOME/.bashrc"
        echo "[ -f \$HOME/.devcontainer_aliases ] && source \$HOME/.devcontainer_aliases" >> "$HOME/.bashrc"
    fi
    
    echo -e "${GREEN}✅ 快捷命令別名已建立${NC}"
    echo -e "${YELLOW}💡 重新啟動 shell 或執行 'source ~/.bashrc' 來啟用別名${NC}"
}

# 檢查並啟動服務
start_services() {
    echo ""
    echo -e "${BLUE}🔌 檢查並啟動服務...${NC}"
    
    # 檢查是否需要啟動 Jupyter
    if [ -f ".devcontainer/.current-profile" ]; then
        local current=$(cat .devcontainer/.current-profile)
        local profile_file=".devcontainer/profiles/${current}.env"
        if [ -f "$profile_file" ]; then
            source "$profile_file"
            if [ "${ENABLE_JUPYTER:-false}" == "true" ]; then
                if ! pgrep -f "jupyter-lab" > /dev/null; then
                    echo "啟動 Jupyter Lab..."
                    nohup jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root > /tmp/jupyter.log 2>&1 &
                    sleep 2
                    echo -e "${GREEN}✅ Jupyter Lab 已啟動 (port 8888)${NC}"
                else
                    echo -e "${GREEN}✅ Jupyter Lab 已在運行${NC}"
                fi
            fi
        fi
    fi
}

# 顯示下一步建議
show_next_steps() {
    echo ""
    echo -e "${CYAN}🎯 建議的下一步操作${NC}"
    echo "========================"
    echo ""
    echo "1. 📝 建立或開啟您的專案檔案"
    echo "2. 🧪 測試環境設定："
    echo "   python -c \"import pandas; print('Python 環境正常!')\""
    echo "   R -e \"library(ggplot2); cat('R 環境正常!\\n')\""
    echo ""
    echo "3. 📚 如果使用 Jupyter，請在瀏覽器開啟："
    echo "   http://localhost:8888"
    echo ""
    echo "4. 🔧 需要額外套件時執行："
    echo "   add-package python <套件名稱>"
    echo "   add-package r <套件名稱>"
    echo ""
}

# 主程式
main() {
    # 檢查是否在正確的目錄
    if [ ! -d ".devcontainer" ]; then
        echo -e "${RED}❌ 錯誤: 請在包含 .devcontainer 目錄的專案根目錄執行此腳本${NC}"
        exit 1
    fi
    
    # 顯示歡迎訊息
    show_welcome
    
    # 檢查現有環境
    if ! check_existing_environment; then
        # 執行環境設定
        setup_environment
    fi
    
    # 顯示環境狀態
    show_status
    
    # 建立快捷命令
    create_aliases
    
    # 啟動服務
    start_services
    
    # 顯示可用工具
    show_available_tools
    
    # 顯示下一步建議
    show_next_steps
    
    echo ""
    echo -e "${GREEN}🎉 Dev Container 快速啟動完成！${NC}"
    echo -e "${CYAN}💡 使用 'status' 命令隨時檢查環境狀態${NC}"
    echo ""
}

# 執行主程式
main "$@"