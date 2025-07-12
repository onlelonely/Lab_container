# 環境備份還原指南

**版本**: v2.0  
**更新日期**: 2025年7月12日  

## 📋 概覽

本指南說明如何備份和還原 Lab Container 環境狀態，包括資料、配置、容器映像檔以及分析結果，確保工作環境的持續性和可復現性。

## 🎯 備份策略

### 完整備份內容
- **環境配置**: Profiles、環境變數、Docker 配置
- **資料檔案**: 原始資料、處理結果、中間檔案
- **容器映像檔**: 自訂映像檔、工具容器
- **腳本和程式碼**: 分析腳本、自訂工具

## 💾 自動備份系統

### 備份腳本設置

```bash
#!/bin/bash
# scripts/backup.sh

BACKUP_DIR="/backups/lab_container"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BACKUP_PATH="${BACKUP_DIR}/${TIMESTAMP}"

echo "開始備份 Lab Container 環境..."
echo "備份路徑: ${BACKUP_PATH}"

# 建立備份目錄
mkdir -p "${BACKUP_PATH}"

# 備份配置檔案
echo "備份配置檔案..."
tar -czf "${BACKUP_PATH}/configs.tar.gz" \
    profiles/ \
    .env \
    scripts/

# 備份資料
echo "備份資料檔案..."
tar -czf "${BACKUP_PATH}/data.tar.gz" \
    --exclude="*.tmp" \
    --exclude="cache/*" \
    data/

echo "備份完成: ${BACKUP_PATH}"
```

### 自動備份排程

```bash
# 設置 crontab 自動備份
crontab -e

# 每日凌晨 2 點進行備份
0 2 * * * /home/cc1296/New_container/Lab_container/scripts/backup.sh
```

## 🔄 還原程序

### 完整環境還原

```bash
#!/bin/bash
# scripts/restore.sh

BACKUP_PATH="$1"

if [ ! -d "$BACKUP_PATH" ]; then
    echo "錯誤: 備份目錄不存在"
    exit 1
fi

echo "開始還原 Lab Container 環境..."

# 還原配置檔案
if [ -f "${BACKUP_PATH}/configs.tar.gz" ]; then
    tar -xzf "${BACKUP_PATH}/configs.tar.gz"
    echo "✅ 配置檔案還原完成"
fi

# 還原資料檔案
if [ -f "${BACKUP_PATH}/data.tar.gz" ]; then
    tar -xzf "${BACKUP_PATH}/data.tar.gz"
    echo "✅ 資料檔案還原完成"
fi

echo "還原完成"
```

---

**相關文檔**:
- [環境設置指南](./environment-setup.md)
- [系統狀態監控](./monitoring.md)