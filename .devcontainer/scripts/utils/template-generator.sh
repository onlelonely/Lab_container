#!/bin/bash
# template-generator.sh - 配置範本生成器

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

# 資訊訊息
info_msg() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

# 生成 devcontainer.json
generate_devcontainer_json() {
    local template_file=".devcontainer/templates/devcontainer.json.template"
    local output_file=".devcontainer/devcontainer.json"
    
    info_msg "生成 devcontainer.json..."
    
    # 設定 R 相關擴展
    local r_extensions=""
    if [ "${ENABLE_RSTUDIO:-false}" == "true" ] || [ "${ENABLE_VSCODE_R:-false}" == "true" ]; then
        r_extensions=',
        "REditorSupport.r",
        "RDebugger.r-debugger"'
    fi
    
    # 設定生物資訊相關擴展
    local bio_extensions=""
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        bio_extensions=',
        "ms-vscode.hexeditor",
        "redhat.vscode-yaml"'
    fi
    
    # 設定 R 相關設定
    local r_settings=""
    if [ "${ENABLE_RSTUDIO:-false}" == "true" ] || [ "${ENABLE_VSCODE_R:-false}" == "true" ]; then
        r_settings=',
        "r.rterm.linux": "/usr/bin/R",
        "r.rterm.option": ["--no-save", "--no-restore"],
        "r.sessionWatcher": true'
    fi
    
    # 替換範本變數
    export R_EXTENSIONS="$r_extensions"
    export BIOINFORMATICS_EXTENSIONS="$bio_extensions"
    export R_SETTINGS="$r_settings"
    
    envsubst < "$template_file" > "$output_file"
    
    success_msg "已生成 devcontainer.json"
}

# 生成 docker-compose.yml
generate_docker_compose() {
    local template_file=".devcontainer/templates/docker-compose.yml.template"
    local output_file=".devcontainer/docker-compose.yml"
    
    info_msg "生成 docker-compose.yml..."
    
    # 設定 RStudio 埠號
    local rstudio_port=""
    if [ "${ENABLE_RSTUDIO:-false}" == "true" ]; then
        rstudio_port='
      - "8787:8787"  # RStudio Server'
    fi
    
    # 設定 GATK 服務
    local gatk_service=""
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        gatk_service='
  gatk-tools:
    image: broadinstitute/gatk:latest
    volumes:
      - ../..:/workspace:cached
      - ./data:/workspace/data
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["bioinformatics", "gatk"]
    environment:
      - JAVA_OPTS=-Xmx8g
    networks:
      - devnet'
    fi
    
    # 設定 ColabFold 服務
    local colabfold_service=""
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        colabfold_service='
  colabfold:
    image: colabfold/colabfold:latest
    volumes:
      - ../..:/workspace:cached
      - ./data:/workspace/data
      - ./models:/workspace/models
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["structure", "colabfold"]
    environment:
      - CUDA_VISIBLE_DEVICES=0
    networks:
      - devnet'
    fi
    
    # 設定 AutoDock 服務
    local autodock_service=""
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        autodock_service='
  autodock:
    image: ccsb/autodock-vina:latest
    volumes:
      - ../..:/workspace:cached
      - ./data:/workspace/data
      - ./results:/workspace/results
    working_dir: /workspace
    profiles: ["molecular", "autodock"]
    networks:
      - devnet'
    fi
    
    # 替換範本變數
    export RSTUDIO_PORT="$rstudio_port"
    export GATK_SERVICE="$gatk_service"
    export COLABFOLD_SERVICE="$colabfold_service"
    export AUTODOCK_SERVICE="$autodock_service"
    
    envsubst < "$template_file" > "$output_file"
    
    success_msg "已生成 docker-compose.yml"
}

# 生成 Dockerfile
generate_dockerfile() {
    local output_file=".devcontainer/Dockerfile.generated"
    
    info_msg "生成 Dockerfile..."
    
    cat > "$output_file" << EOF
# Auto-generated Dockerfile for profile: ${PROFILE_NAME}
# Generated on: $(date)

FROM mcr.microsoft.com/devcontainers/anaconda:0-3

# 環境變數
ENV PROFILE=${PROFILE_NAME}
ENV ENABLE_BIOINFORMATICS=${ENABLE_BIOINFORMATICS:-false}

# 安裝系統套件
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \\
    && apt-get -y install --no-install-recommends \\
        build-essential \\
        curl \\
        git \\
        vim \\
        wget \\
        unzip \\
        && apt-get clean -y && rm -rf /var/lib/apt/lists/*

EOF

    # 如果啟用生物資訊工具，添加相關安裝
    if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
        cat >> "$output_file" << EOF
# 安裝生物資訊工具
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \\
    && apt-get -y install --no-install-recommends \\
        default-jre \\
        && apt-get clean -y && rm -rf /var/lib/apt/lists/*

EOF
    fi

    # 如果啟用 RStudio，添加安裝
    if [ "${ENABLE_RSTUDIO:-false}" == "true" ]; then
        cat >> "$output_file" << EOF
# 安裝 RStudio Server
RUN wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2023.12.1-402-amd64.deb \\
    && dpkg -i rstudio-server-2023.12.1-402-amd64.deb || true \\
    && apt-get -f install -y \\
    && rm rstudio-server-2023.12.1-402-amd64.deb

EOF
    fi

    cat >> "$output_file" << EOF
# 複製環境配置
COPY .devcontainer/configs/ /tmp/configs/

# 安裝 conda 環境
RUN if [ -f "/tmp/configs/conda/${CONDA_ENV_FILE}" ]; then \\
        conda env update --file /tmp/configs/conda/${CONDA_ENV_FILE} --name base; \\
    fi

# 安裝 Python 套件
RUN if [ -f "/tmp/configs/python/${PYTHON_REQUIREMENTS_FILE}" ]; then \\
        pip install -r /tmp/configs/python/${PYTHON_REQUIREMENTS_FILE}; \\
    fi

# 安裝 R 套件
RUN if [ -f "/tmp/configs/r/${R_SCRIPT_FILE}" ]; then \\
        Rscript /tmp/configs/r/${R_SCRIPT_FILE}; \\
    fi

# 設定使用者
USER vscode

# 設定工作目錄
WORKDIR /workspace

EOF

    success_msg "已生成 Dockerfile"
}

# 生成 post-create 腳本
generate_post_create() {
    local output_file=".devcontainer/scripts/utils/post-create.sh"
    
    info_msg "生成 post-create 腳本..."
    
    cat > "$output_file" << 'EOF'
#!/bin/bash
# post-create.sh - 容器建立後執行的腳本

set -euo pipefail

echo "🚀 執行 post-create 設定..."

# 設定別名
if [ -f "$HOME/.devcontainer_aliases" ]; then
    source "$HOME/.devcontainer_aliases"
fi

# 檢查環境
echo "📊 環境檢查:"
echo "Profile: ${PROFILE:-未設定}"
echo "Python: $(python --version 2>/dev/null || echo '未安裝')"
echo "R: $(R --version 2>/dev/null | head -1 || echo '未安裝')"

# 如果啟用 Jupyter，設定 kernel
if [ "${ENABLE_JUPYTER:-false}" == "true" ]; then
    echo "🔧 設定 Jupyter kernel..."
    python -m ipykernel install --user --name="python3" --display-name="Python 3"
fi

# 如果啟用生物資訊工具，檢查工具
if [ "${ENABLE_BIOINFORMATICS:-false}" == "true" ]; then
    echo "🧬 檢查生物資訊工具:"
    for tool in fastqc samtools bcftools bedtools; do
        if command -v "$tool" >/dev/null 2>&1; then
            echo "  ✅ $tool"
        else
            echo "  ❌ $tool (未安裝)"
        fi
    done
fi

echo "✅ post-create 設定完成"
EOF
    
    chmod +x "$output_file"
    success_msg "已生成 post-create 腳本"
}

# 主要生成函數
generate_all_configs() {
    if [ ! -f ".devcontainer/.current-profile" ]; then
        error_exit "找不到當前 profile，請先執行環境設定"
    fi
    
    local current_profile=$(cat .devcontainer/.current-profile)
    local profile_file=".devcontainer/profiles/${current_profile}.env"
    
    if [ ! -f "$profile_file" ]; then
        error_exit "找不到 profile 檔案: $profile_file"
    fi
    
    # 載入 profile 設定
    source "$profile_file"
    
    info_msg "為 profile '$PROFILE_NAME' 生成配置檔案..."
    
    # 生成各種配置檔案
    generate_devcontainer_json
    generate_docker_compose
    generate_dockerfile
    generate_post_create
    
    success_msg "所有配置檔案生成完成！"
    echo ""
    echo "生成的檔案："
    echo "  📄 .devcontainer/devcontainer.json"
    echo "  📄 .devcontainer/docker-compose.yml"
    echo "  📄 .devcontainer/Dockerfile.generated"
    echo "  📄 .devcontainer/scripts/utils/post-create.sh"
}

# 主程式
main() {
    # 檢查是否在正確的目錄
    if [ ! -d ".devcontainer" ]; then
        error_exit "請在包含 .devcontainer 目錄的專案根目錄執行此腳本"
    fi
    
    generate_all_configs
}

# 執行主程式
main "$@"