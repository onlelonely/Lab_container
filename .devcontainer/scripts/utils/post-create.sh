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
