# 系統監控指南

**文檔版本**: v2.0  
**適用系統**: Lab Container v2.0+  
**目標讀者**: 系統管理員、進階使用者

## 📊 監控概述

Lab Container 提供多層次的監控系統，幫助您了解系統狀態、效能表現和資源使用情況。

### 監控層級
```
監控系統架構:
├── 🖥️ 系統層監控 (主機資源)
├── 🐳 容器層監控 (Docker 容器)
├── 📦 應用層監控 (工具與服務)
└── 📈 業務層監控 (分析工作流程)
```

## 🛠 內建監控工具

### 系統狀態檢查器
```bash
# 完整系統狀態檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh

# 詳細模式
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose

# 指定組件檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh --component docker
bash .devcontainer/scripts/manage/devcontainer-status.sh --component tools
bash .devcontainer/scripts/manage/devcontainer-status.sh --component packages
```

### 工具可用性監控
```bash
# 檢查所有工具狀態
bash .devcontainer/scripts/manage/tool-manager.sh check all

# 檢查特定工具
bash .devcontainer/scripts/manage/tool-manager.sh check gatk
bash .devcontainer/scripts/manage/tool-manager.sh check colabfold
bash .devcontainer/scripts/manage/tool-manager.sh check autodock

# 工具健康檢查
bash .devcontainer/scripts/manage/tool-manager.sh health
```

## 📈 資源監控

### 即時資源監控
```bash
# Docker 容器資源使用
docker stats

# 指定容器監控
docker stats lab-container-devcontainer-1

# 無流式輸出 (一次性)
docker stats --no-stream

# 格式化輸出
docker stats --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.NetIO}}\t{{.BlockIO}}"
```

### 磁碟使用監控
```bash
# 系統磁碟使用
df -h

# Docker 系統使用
docker system df

# 詳細 Docker 使用資訊
docker system df -v

# 專案目錄使用
du -sh data/ results/ models/ temp/ backups/

# 深入分析大檔案
find . -type f -size +100M -exec ls -lh {} \; | sort -k5 -hr
```

### 記憶體使用監控
```bash
# 系統記憶體概況
free -h

# 詳細記憶體資訊
cat /proc/meminfo

# 進程記憶體使用
ps aux --sort=-%mem | head -10

# 容器記憶體限制檢查
docker inspect lab-container-devcontainer-1 | grep -i memory
```

## 🔍 日誌監控

### 容器日誌收集
```bash
# 主開發容器日誌
docker logs lab-container-devcontainer-1

# 即時日誌監控
docker logs -f lab-container-devcontainer-1

# 最近 N 行日誌
docker logs --tail 100 lab-container-devcontainer-1

# 時間範圍日誌
docker logs --since="2025-07-12T10:00:00" lab-container-devcontainer-1
```

### 工具容器日誌
```bash
# 查看 GATK 容器日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk

# 即時監控工具日誌
bash .devcontainer/scripts/manage/tool-manager.sh manage logs gatk -f

# ColabFold 日誌 (如果正在運行)
docker logs lab-container-colabfold-1

# 所有工具容器日誌
for container in $(docker ps -a --format "{{.Names}}" | grep lab-container); do
    echo "=== $container ==="
    docker logs --tail 20 $container
done
```

### 應用程式日誌
```bash
# Jupyter Lab 日誌
tail -f ~/.local/share/jupyter/runtime/jupyter.log

# R 會話日誌
tail -f ~/.RData

# Conda 日誌
tail -f ~/.conda/envs/*/conda-meta/*.json
```

## 📊 效能監控

### CPU 效能監控
```bash
# 即時 CPU 使用率
top -p $(docker inspect --format='{{.State.Pid}}' lab-container-devcontainer-1)

# htop 互動式監控 (如已安裝)
htop -p $(docker inspect --format='{{.State.Pid}}' lab-container-devcontainer-1)

# CPU 使用歷史
sar -u 1 5  # 每秒採樣，共5次

# 負載平均
uptime
```

### I/O 效能監控
```bash
# 磁碟 I/O 統計
iostat -x 1

# 進程 I/O 監控
iotop -o

# 檔案系統 I/O
sar -d 1 5

# 容器 I/O 統計
docker stats --format "table {{.Container}}\t{{.BlockIO}}\t{{.NetIO}}"
```

### 網路效能監控
```bash
# 網路介面統計
netstat -i

# 網路連線監控
ss -tuln

# 容器網路統計
docker exec lab-container-devcontainer-1 cat /proc/net/dev

# 網路速度測試
speedtest-cli  # 如已安裝
```

## 🚨 警報與通知

### 資源警報腳本
```bash
#!/bin/bash
# resource-alert.sh - 資源使用警報腳本

# 設定閾值
CPU_THRESHOLD=80
MEMORY_THRESHOLD=80
DISK_THRESHOLD=90

# 檢查 CPU 使用率
CPU_USAGE=$(docker stats --no-stream --format "{{.CPUPerc}}" lab-container-devcontainer-1 | sed 's/%//')
if (( $(echo "$CPU_USAGE > $CPU_THRESHOLD" | bc -l) )); then
    echo "⚠️  CPU 使用率過高: ${CPU_USAGE}%"
    # 可在此加入通知邏輯 (email, slack 等)
fi

# 檢查記憶體使用率
MEMORY_USAGE=$(docker stats --no-stream --format "{{.MemPerc}}" lab-container-devcontainer-1 | sed 's/%//')
if (( $(echo "$MEMORY_USAGE > $MEMORY_THRESHOLD" | bc -l) )); then
    echo "⚠️  記憶體使用率過高: ${MEMORY_USAGE}%"
fi

# 檢查磁碟使用率
DISK_USAGE=$(df / | awk 'NR==2 {print $5}' | sed 's/%//')
if [ "$DISK_USAGE" -gt "$DISK_THRESHOLD" ]; then
    echo "⚠️  磁碟使用率過高: ${DISK_USAGE}%"
fi
```

### 自動監控腳本
```bash
#!/bin/bash
# auto-monitor.sh - 自動監控腳本

LOG_FILE="monitoring/$(date +%Y%m%d)_monitor.log"
mkdir -p monitoring

while true; do
    TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
    
    # 系統資源
    CPU_USAGE=$(docker stats --no-stream --format "{{.CPUPerc}}" lab-container-devcontainer-1)
    MEMORY_USAGE=$(docker stats --no-stream --format "{{.MemUsage}}" lab-container-devcontainer-1)
    
    # 記錄日誌
    echo "[$TIMESTAMP] CPU: $CPU_USAGE, Memory: $MEMORY_USAGE" >> "$LOG_FILE"
    
    # 檢查服務狀態
    if ! docker ps | grep -q lab-container-devcontainer-1; then
        echo "[$TIMESTAMP] ❌ 主容器未運行" >> "$LOG_FILE"
    fi
    
    # 每30秒檢查一次
    sleep 30
done
```

## 📈 效能基準測試

### 容器啟動時間測試
```bash
#!/bin/bash
# startup-benchmark.sh - 容器啟動效能測試

PROFILES=("minimal" "datascience" "bioinformatics" "ml" "statistics" "full")

for profile in "${PROFILES[@]}"; do
    echo "🧪 測試 $profile Profile 啟動時間..."
    
    start_time=$(date +%s)
    
    # 設置環境
    bash .devcontainer/scripts/manage/setup-environment.sh "$profile" > /dev/null 2>&1
    
    # 重建容器 (模擬)
    docker-compose -f .devcontainer/docker-compose.yml down > /dev/null 2>&1
    docker-compose -f .devcontainer/docker-compose.yml up -d > /dev/null 2>&1
    
    # 等待容器就緒
    until docker exec lab-container-devcontainer-1 echo "ready" > /dev/null 2>&1; do
        sleep 1
    done
    
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    echo "✅ $profile Profile 啟動耗時: ${duration} 秒"
done
```

### 套件安裝效能測試
```bash
#!/bin/bash
# package-install-benchmark.sh - 套件安裝效能測試

TEST_PACKAGES=("pandas" "numpy" "matplotlib" "seaborn" "scikit-learn")

for package in "${TEST_PACKAGES[@]}"; do
    echo "🧪 測試 $package 安裝時間..."
    
    start_time=$(date +%s)
    bash .devcontainer/scripts/utils/package-manager.sh add python "$package" > /dev/null 2>&1
    end_time=$(date +%s)
    
    duration=$((end_time - start_time))
    echo "✅ $package 安裝耗時: ${duration} 秒"
done
```

## 📋 監控儀表板

### 簡易文字儀表板
```bash
#!/bin/bash
# dashboard.sh - 簡易監控儀表板

clear
echo "======================================"
echo "    Lab Container 監控儀表板"
echo "======================================"
echo

# 系統資訊
echo "🖥️  系統資訊:"
echo "   時間: $(date)"
echo "   主機: $(hostname)"
echo "   運行時間: $(uptime | awk '{print $3,$4}' | sed 's/,//')"
echo

# 容器狀態
echo "🐳 容器狀態:"
docker ps --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" | grep lab-container
echo

# 資源使用
echo "📊 資源使用:"
docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"
echo

# 磁碟使用
echo "💾 磁碟使用:"
df -h | grep -E "(Filesystem|/dev/)"
echo

# 工具狀態
echo "🧬 工具狀態:"
for tool in gatk colabfold autodock; do
    if docker ps | grep -q "$tool"; then
        echo "   $tool: ✅ 運行中"
    else
        echo "   $tool: ⏸️  停止"
    fi
done
```

### Web 儀表板 (使用 Grafana + Prometheus)
```yaml
# docker-compose.monitoring.yml - 監控服務
version: '3.8'

services:
  prometheus:
    image: prom/prometheus:latest
    ports:
      - "9090:9090"
    volumes:
      - ./monitoring/prometheus.yml:/etc/prometheus/prometheus.yml
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.path=/prometheus'
      - '--web.console.libraries=/etc/prometheus/console_libraries'
      - '--web.console.templates=/etc/prometheus/consoles'

  grafana:
    image: grafana/grafana:latest
    ports:
      - "3000:3000"
    environment:
      - GF_SECURITY_ADMIN_PASSWORD=admin
    volumes:
      - grafana_data:/var/lib/grafana

  cadvisor:
    image: gcr.io/cadvisor/cadvisor:latest
    ports:
      - "8080:8080"
    volumes:
      - /:/rootfs:ro
      - /var/run:/var/run:ro
      - /sys:/sys:ro
      - /var/lib/docker/:/var/lib/docker:ro

volumes:
  grafana_data:
```

## 🎯 監控最佳實踐

### 1. 定期監控
```bash
# 設定 crontab 定期監控
# 每5分鐘檢查系統狀態
*/5 * * * * /path/to/resource-alert.sh

# 每小時生成效能報告
0 * * * * /path/to/performance-report.sh

# 每天清理舊日誌
0 2 * * * find /path/to/logs -name "*.log" -mtime +7 -delete
```

### 2. 閾值設定
```bash
# 建議的監控閾值
CPU_WARNING=70        # CPU 使用率警告閾值
CPU_CRITICAL=90       # CPU 使用率臨界閾值

MEMORY_WARNING=80     # 記憶體使用率警告閾值
MEMORY_CRITICAL=95    # 記憶體使用率臨界閾值

DISK_WARNING=85       # 磁碟使用率警告閾值
DISK_CRITICAL=95      # 磁碟使用率臨界閾值

LOAD_WARNING=2.0      # 負載平均警告閾值
LOAD_CRITICAL=5.0     # 負載平均臨界閾值
```

### 3. 日誌輪替
```bash
# 設定日誌輪替 /etc/logrotate.d/labcontainer
/path/to/lab-container/logs/*.log {
    daily
    missingok
    rotate 30
    compress
    notifempty
    create 644 user user
}
```

### 4. 自動恢復
```bash
#!/bin/bash
# auto-recovery.sh - 自動恢復腳本

# 檢查主容器狀態
if ! docker ps | grep -q lab-container-devcontainer-1; then
    echo "⚠️  主容器異常，嘗試重啟..."
    docker-compose -f .devcontainer/docker-compose.yml restart devcontainer
fi

# 檢查磁碟空間
DISK_USAGE=$(df / | awk 'NR==2 {print $5}' | sed 's/%//')
if [ "$DISK_USAGE" -gt 90 ]; then
    echo "⚠️  磁碟空間不足，執行清理..."
    docker system prune -f
fi

# 檢查記憶體使用
MEMORY_USAGE=$(free | awk 'NR==2{printf "%.2f", $3*100/$2}')
if (( $(echo "$MEMORY_USAGE > 90" | bc -l) )); then
    echo "⚠️  記憶體使用過高，重啟服務..."
    docker-compose -f .devcontainer/docker-compose.yml restart
fi
```

## 📊 報告生成

### 每日效能報告
```bash
#!/bin/bash
# daily-report.sh - 每日效能報告

DATE=$(date +%Y-%m-%d)
REPORT_FILE="reports/daily_report_$DATE.txt"
mkdir -p reports

{
    echo "Lab Container 每日效能報告"
    echo "日期: $DATE"
    echo "=================================="
    echo
    
    echo "📊 系統資源使用:"
    docker stats --no-stream
    echo
    
    echo "💾 磁碟使用情況:"
    df -h
    echo
    
    echo "🐳 容器狀態:"
    docker ps -a
    echo
    
    echo "🧬 工具可用性:"
    bash .devcontainer/scripts/manage/tool-manager.sh check all
    echo
    
    echo "📈 當日活動總結:"
    echo "- 容器重啟次數: $(docker events --since="$DATE 00:00:00" --until="$DATE 23:59:59" --filter event=restart | wc -l)"
    echo "- 錯誤日誌條目: $(docker logs lab-container-devcontainer-1 --since="$DATE 00:00:00" 2>&1 | grep -i error | wc -l)"
    
} > "$REPORT_FILE"

echo "📋 每日報告已生成: $REPORT_FILE"
```

---

**有效的監控是系統穩定運行的基礎**  
*持續監控，預防問題，確保最佳效能*