#!/bin/bash
# tool-manager.sh - 專用工具容器管理系統

set -euo pipefail

# 顏色定義
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
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

# 檢查 Docker 和 Docker Compose
check_docker() {
    if ! command -v docker >/dev/null 2>&1; then
        error_exit "Docker 未安裝或不可用"
    fi
    
    if ! command -v docker-compose >/dev/null 2>&1; then
        error_exit "Docker Compose 未安裝或不可用"
    fi
    
    if ! docker ps >/dev/null 2>&1; then
        error_exit "Docker daemon 未運行或權限不足"
    fi
}

# 檢查 GPU 可用性
check_gpu() {
    if command -v nvidia-smi >/dev/null 2>&1; then
        if nvidia-smi >/dev/null 2>&1; then
            return 0  # GPU 可用
        fi
    fi
    return 1  # GPU 不可用
}

# GATK 工具包裝器
run_gatk() {
    info_msg "執行 GATK: $*"
    
    if ! docker-compose ps gatk-tools >/dev/null 2>&1; then
        info_msg "啟動 GATK 容器..."
        docker-compose --profile gatk up -d gatk-tools
        sleep 5
    fi
    
    docker-compose --profile gatk exec gatk-tools gatk "$@"
}

# ColabFold 工具包裝器
run_colabfold() {
    info_msg "執行 ColabFold: $*"
    
    # 檢查 GPU 可用性
    if ! check_gpu; then
        warning_msg "未檢測到 GPU，將使用 CPU 模式（速度較慢）"
        docker-compose run --rm colabfold colabfold_batch --num-models 1 "$@"
    else
        success_msg "檢測到 GPU，使用 GPU 加速模式"
        docker-compose --profile colabfold run --rm colabfold colabfold_batch "$@"
    fi
}

# AutoDock Vina 工具包裝器
run_autodock() {
    info_msg "執行 AutoDock Vina: $*"
    
    if ! docker-compose ps autodock >/dev/null 2>&1; then
        info_msg "啟動 AutoDock 容器..."
        docker-compose --profile autodock up -d autodock
        sleep 3
    fi
    
    docker-compose --profile autodock exec autodock vina "$@"
}

# 工具狀態檢查
check_tool_availability() {
    local tool="$1"
    
    echo -e "${CYAN}🔍 檢查工具: $tool${NC}"
    echo "======================"
    
    case $tool in
        "gatk")
            if docker-compose --profile gatk run --rm gatk-tools gatk --version 2>/dev/null; then
                success_msg "GATK 可用"
            else
                error_exit "GATK 不可用"
            fi
            ;;
        "colabfold")
            if docker-compose --profile colabfold run --rm colabfold colabfold_batch --help >/dev/null 2>&1; then
                success_msg "ColabFold 可用"
                if check_gpu; then
                    success_msg "GPU 加速可用"
                else
                    warning_msg "僅 CPU 模式可用"
                fi
            else
                error_exit "ColabFold 不可用"
            fi
            ;;
        "autodock")
            if docker-compose --profile autodock run --rm autodock vina --version 2>/dev/null; then
                success_msg "AutoDock Vina 可用"
            else
                error_exit "AutoDock Vina 不可用"
            fi
            ;;
        "all")
            check_tool_availability "gatk"
            echo ""
            check_tool_availability "colabfold"
            echo ""
            check_tool_availability "autodock"
            ;;
        *)
            error_exit "未知工具: $tool。可用工具: gatk, colabfold, autodock, all"
            ;;
    esac
}

# 批次工具管理
manage_tools() {
    local action="$1"
    local tools="${2:-all}"
    
    case $action in
        "start")
            info_msg "啟動工具容器..."
            case $tools in
                "gatk")
                    docker-compose --profile gatk up -d gatk-tools
                    ;;
                "colabfold")
                    docker-compose --profile colabfold up -d colabfold
                    ;;
                "autodock")
                    docker-compose --profile autodock up -d autodock
                    ;;
                "all"|*)
                    docker-compose --profile gatk --profile colabfold --profile autodock up -d
                    ;;
            esac
            success_msg "工具容器啟動完成"
            ;;
        "stop")
            info_msg "停止工具容器..."
            case $tools in
                "gatk")
                    docker-compose --profile gatk down
                    ;;
                "colabfold")
                    docker-compose --profile colabfold down
                    ;;
                "autodock")
                    docker-compose --profile autodock down
                    ;;
                "all"|*)
                    docker-compose --profile gatk --profile colabfold --profile autodock down
                    ;;
            esac
            success_msg "工具容器停止完成"
            ;;
        "restart")
            manage_tools "stop" "$tools"
            sleep 2
            manage_tools "start" "$tools"
            ;;
        "status")
            echo -e "${CYAN}📊 工具容器狀態${NC}"
            echo "===================="
            echo ""
            docker-compose ps
            
            echo ""
            echo -e "${CYAN}💾 容器資源使用${NC}"
            echo "===================="
            docker stats --no-stream --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}" | grep -E "(gatk|colabfold|autodock|NAME)" || echo "無相關容器運行"
            ;;
        "logs")
            local service="${tools:-gatk-tools}"
            info_msg "顯示 $service 日誌..."
            docker-compose logs --tail=50 "$service"
            ;;
        "clean")
            warning_msg "清理未使用的容器和映像..."
            docker system prune -f
            docker volume prune -f
            success_msg "清理完成"
            ;;
        *)
            error_exit "未知動作: $action。可用動作: start, stop, restart, status, logs, clean"
            ;;
    esac
}

# 建立工作目錄
setup_directories() {
    info_msg "建立工作目錄..."
    mkdir -p data results models
    
    # 建立範例資料目錄結構
    mkdir -p data/{fastq,reference,variants}
    mkdir -p results/{qc,alignment,variants,structure,docking}
    mkdir -p models/{colabfold,autodock}
    
    success_msg "工作目錄建立完成"
    echo "目錄結構："
    echo "  📁 data/ - 輸入資料"
    echo "    📁 fastq/ - FASTQ 檔案"
    echo "    📁 reference/ - 參考基因組"
    echo "    📁 variants/ - 變異檔案"
    echo "  📁 results/ - 分析結果"
    echo "    📁 qc/ - 品質控制結果"
    echo "    📁 alignment/ - 比對結果"
    echo "    📁 variants/ - 變異分析結果"
    echo "    📁 structure/ - 結構預測結果"
    echo "    📁 docking/ - 分子對接結果"
    echo "  📁 models/ - 模型檔案"
}

# 工作流程範例
show_workflow_examples() {
    echo -e "${CYAN}🧬 工作流程範例${NC}"
    echo "=================="
    echo ""
    
    echo -e "${GREEN}📊 NGS 變異分析工作流程${NC}"
    echo "1. 品質控制: fastqc data/fastq/*.fastq.gz -o results/qc/"
    echo "2. 序列比對: bwa mem + samtools"
    echo "3. 變異呼叫: $0 gatk HaplotypeCaller -I sample.bam -O results/variants/variants.vcf -R reference.fa"
    echo "4. 變異過濾: $0 gatk VariantFiltration"
    echo ""
    
    echo -e "${GREEN}🧬 蛋白質結構預測工作流程${NC}"
    echo "1. 準備序列: 將蛋白質序列儲存為 FASTA 格式"
    echo "2. 結構預測: $0 colabfold protein.fasta results/structure/ --num-models 5"
    echo "3. 結果分析: 檢視 PDB 檔案和信心分數"
    echo ""
    
    echo -e "${GREEN}💊 分子對接工作流程${NC}"
    echo "1. 準備受體: 蛋白質結構 PDB 檔案"
    echo "2. 準備配體: 小分子 SDF 或 MOL2 檔案"
    echo "3. 分子對接: $0 autodock --receptor protein.pdbqt --ligand ligand.pdbqt --out results/docking/"
    echo ""
}

# 顯示使用說明
show_help() {
    echo ""
    echo -e "${BLUE}🧬 生物資訊工具容器管理器${NC}"
    echo "============================="
    echo ""
    echo -e "${GREEN}輕量級工具 (dev container 內)${NC}"
    echo "  fastqc <files>     - 品質控制分析"
    echo "  samtools <cmd>     - SAM/BAM 操作"
    echo "  bcftools <cmd>     - VCF 操作"
    echo "  bedtools <cmd>     - 基因組區間操作"
    echo ""
    echo -e "${GREEN}重量級工具 (專用容器)${NC}"
    echo "  $0 gatk <args>     - GATK 變異呼叫"
    echo "  $0 colabfold <args> - 蛋白質結構預測"
    echo "  $0 autodock <args> - 分子對接"
    echo ""
    echo -e "${GREEN}管理命令${NC}"
    echo "  $0 check <tool>    - 檢查工具可用性"
    echo "  $0 manage <action> [tool] - 管理工具容器"
    echo "  $0 setup           - 建立工作目錄"
    echo "  $0 workflow        - 顯示工作流程範例"
    echo ""
    echo -e "${GREEN}管理動作${NC}"
    echo "  start    - 啟動容器"
    echo "  stop     - 停止容器"
    echo "  restart  - 重啟容器"
    echo "  status   - 檢查狀態"
    echo "  logs     - 顯示日誌"
    echo "  clean    - 清理資源"
    echo ""
    echo -e "${GREEN}範例${NC}"
    echo "  $0 gatk --version"
    echo "  $0 colabfold protein.fasta results/structure/"
    echo "  $0 check all"
    echo "  $0 manage start gatk"
    echo "  $0 manage status"
    echo ""
}

# 主程式
main() {
    # 檢查 Docker 環境
    check_docker
    
    # 檢查參數
    if [ $# -eq 0 ]; then
        show_help
        return 0
    fi
    
    local command="$1"
    shift
    
    case "$command" in
        "gatk")
            run_gatk "$@"
            ;;
        "colabfold")
            run_colabfold "$@"
            ;;
        "autodock")
            run_autodock "$@"
            ;;
        "check")
            local tool="${1:-all}"
            check_tool_availability "$tool"
            ;;
        "manage")
            if [ $# -lt 1 ]; then
                error_exit "用法: $0 manage <action> [tool]"
            fi
            local action="$1"
            local tool="${2:-all}"
            manage_tools "$action" "$tool"
            ;;
        "setup")
            setup_directories
            ;;
        "workflow")
            show_workflow_examples
            ;;
        "help"|"-h"|"--help")
            show_help
            ;;
        *)
            error_exit "未知命令: $command。使用 '$0 help' 查看使用說明"
            ;;
    esac
}

# 執行主程式
main "$@"