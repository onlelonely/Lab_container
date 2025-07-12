# 故障排除指南

**文檔版本**: v2.0  
**適用系統**: Lab Container v2.0+  
**目標讀者**: 所有使用者

## 🚨 常見問題與解決方案

### 🐳 Docker 相關問題

#### 問題: Docker 服務未啟動
```bash
# 錯誤訊息
Cannot connect to the Docker daemon at unix:///var/run/docker.sock

# 解決方案
# Ubuntu/Debian
sudo systemctl start docker
sudo systemctl enable docker

# macOS (Docker Desktop)
# 確保 Docker Desktop 應用程式正在運行

# Windows (Docker Desktop)
# 確保 Docker Desktop 應用程式正在運行且 WSL2 已啟用
```

#### 問題: 權限被拒絕
```bash
# 錯誤訊息
Got permission denied while trying to connect to the Docker daemon socket

# 解決方案
# 將使用者加入 docker 群組
sudo usermod -aG docker $USER
# 重新登入或執行
newgrp docker

# 測試
docker run hello-world
```

#### 問題: 磁碟空間不足
```bash
# 錯誤訊息
No space left on device

# 診斷
df -h
docker system df

# 解決方案
# 清理不使用的容器、映像和快取
docker system prune -a
docker volume prune

# 清理特定組件
docker container prune  # 移除停止的容器
docker image prune -a   # 移除未使用的映像
docker network prune    # 移除未使用的網路
```

### 📦 套件安裝問題

#### 問題: Python 套件安裝失敗
```bash
# 錯誤訊息範例
ERROR: Could not build wheels for package_name

# 解決方案 1: 使用 Conda 安裝
bash .devcontainer/scripts/utils/package-manager.sh add conda package_name

# 解決方案 2: 安裝編譯依賴
sudo apt-get update
sudo apt-get install build-essential python3-dev

# 解決方案 3: 使用預編譯版本
pip install --only-binary=all package_name

# 解決方案 4: 增加暫存空間
export TMPDIR=/tmp/pip-build
mkdir -p $TMPDIR
pip install package_name
```

#### 問題: R 套件安裝失敗
```bash
# 錯誤訊息範例
compilation of package 'package_name' failed

# 解決方案 1: 安裝系統依賴
sudo apt-get install r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

# 解決方案 2: 使用 Bioconductor
R -e "BiocManager::install('package_name')"

# 解決方案 3: 安裝二進位套件
R -e "install.packages('package_name', type='binary')"

# 解決方案 4: 從 CRAN 快照安裝
R -e "install.packages('package_name', repos='https://cran.microsoft.com/snapshot/2024-01-01')"
```

#### 問題: Conda 環境損壞
```bash
# 症狀
conda: command not found
或套件版本衝突

# 解決方案 1: 重新載入 Conda 環境
source /opt/conda/bin/activate
conda init bash
source ~/.bashrc

# 解決方案 2: 修復環境
conda update conda
conda clean --all
conda env update --file .devcontainer/configs/conda/environment-full.yml

# 解決方案 3: 重建環境 (最後手段)
conda env remove -n base
bash .devcontainer/scripts/manage/setup-environment.sh <profile_name>
```

### 🔧 容器建構問題

#### 問題: DevContainer 建構失敗
```bash
# 錯誤訊息
Failed to create container

# 診斷步驟
# 1. 檢查配置語法
cat .devcontainer/devcontainer.json | jq .

# 2. 檢查 Docker Compose 語法  
docker-compose -f .devcontainer/docker-compose.yml config

# 3. 手動建構測試
docker-compose -f .devcontainer/docker-compose.yml build

# 4. 檢查日誌
docker-compose -f .devcontainer/docker-compose.yml logs devcontainer
```

#### 問題: 基礎映像下載失敗
```bash
# 錯誤訊息
Error response from daemon: pull access denied

# 解決方案 1: 檢查網路連線
ping hub.docker.com
curl -I https://registry-1.docker.io/

# 解決方案 2: 使用國內映像源
# 編輯 /etc/docker/daemon.json
{
  "registry-mirrors": [
    "https://registry.cn-hangzhou.aliyuncs.com",
    "https://hub-mirror.c.163.com"
  ]
}

# 重啟 Docker
sudo systemctl restart docker

# 解決方案 3: 手動拉取
docker pull mcr.microsoft.com/devcontainers/anaconda:0-3
```

### 🧬 生物資訊工具問題

#### 問題: GATK 容器啟動失敗
```bash
# 診斷
docker logs lab-container-gatk-tools-1

# 常見解決方案
# 1. 記憶體不足
# 編輯 docker-compose.tools.yml 增加記憶體限制
deploy:
  resources:
    limits:
      memory: 16G

# 2. GPU 不可用 (如果使用 GPU 版本)
nvidia-smi
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi

# 3. 重新拉取映像
docker pull broadinstitute/gatk:latest
bash .devcontainer/scripts/manage/tool-manager.sh manage restart gatk
```

#### 問題: ColabFold GPU 不工作
```bash
# 診斷 GPU 狀態
nvidia-smi
nvidia-docker version

# 解決方案 1: 檢查 NVIDIA Container Toolkit
sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker

# 解決方案 2: 確認 Docker GPU 支援
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi

# 解決方案 3: 檢查容器 GPU 配置
# 確保 docker-compose.tools.yml 中包含:
runtime: nvidia
environment:
  - NVIDIA_VISIBLE_DEVICES=all
```

### 🌐 網路連線問題

#### 問題: 無法訪問外部套件庫
```bash
# 症狀
Connection timeout
SSL certificate problem

# 解決方案 1: 檢查 DNS
nslookup pypi.org
nslookup cran.r-project.org

# 解決方案 2: 設定代理 (如需要)
export http_proxy=http://proxy:port
export https_proxy=https://proxy:port

# 解決方案 3: 使用國內源
# Python
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple/ package_name

# R
options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Conda
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
```

#### 問題: 容器間通訊失敗
```bash
# 診斷網路
docker network ls
docker network inspect lab-container_devnet

# 解決方案 1: 重建網路
docker-compose -f .devcontainer/docker-compose.yml down
docker network prune
docker-compose -f .devcontainer/docker-compose.yml up -d

# 解決方案 2: 檢查防火牆設定
sudo ufw status
sudo iptables -L DOCKER-USER
```

### 💾 資料存取問題

#### 問題: Volume 掛載失敗
```bash
# 錯誤訊息
bind: no such file or directory

# 解決方案 1: 檢查路徑存在
ls -la ./data ./results

# 解決方案 2: 建立缺失目錄
mkdir -p data results scripts models temp backups

# 解決方案 3: 檢查權限
chmod 755 data results
chown -R $USER:$USER data results

# 解決方案 4: 絕對路徑掛載
# 編輯 docker-compose.yml 使用絕對路徑
volumes:
  - $(pwd)/data:/workspace/data
```

#### 問題: 檔案權限錯誤
```bash
# 症狀
Permission denied

# 解決方案 1: 檢查檔案擁有者
ls -la problematic_file

# 解決方案 2: 修改權限
chmod 644 data_file
chmod 755 script_file

# 解決方案 3: 變更擁有者
sudo chown vscode:vscode file_or_directory

# 解決方案 4: 容器內修正
docker exec -it container_name chown -R vscode:vscode /workspace
```

### 🔄 環境同步問題

#### 問題: Profile 切換後環境不一致
```bash
# 症狀
套件版本不符預期
工具不可用

# 解決方案 1: 強制重新生成
bash .devcontainer/scripts/manage/setup-environment.sh <profile> --force

# 解決方案 2: 清理快取
rm -rf .devcontainer/.generated/
docker system prune -f

# 解決方案 3: 檢查 Profile 配置
cat .devcontainer/profiles/<profile>.env
bash -n .devcontainer/profiles/<profile>.env
```

#### 問題: 套件版本衝突
```bash
# 診斷衝突
pip check
conda list --show-channel-urls

# 解決方案 1: 重新安裝衝突套件
pip install --force-reinstall conflicting_package

# 解決方案 2: 使用虛擬環境
python -m venv temp_env
source temp_env/bin/activate
pip install problem_package

# 解決方案 3: 降級相依套件
pip install "dependency_package<version"
```

## 🛠 診斷工具

### 系統狀態檢查
```bash
# 完整系統診斷
bash .devcontainer/scripts/manage/devcontainer-status.sh

# 詳細檢查
bash .devcontainer/scripts/manage/devcontainer-status.sh --verbose

# 檢查特定組件
bash .devcontainer/scripts/manage/tool-manager.sh check all
```

### 日誌收集
```bash
# 收集所有日誌
mkdir -p troubleshooting/logs
docker-compose -f .devcontainer/docker-compose.yml logs > troubleshooting/logs/compose.log
docker logs lab-container-devcontainer-1 > troubleshooting/logs/devcontainer.log

# 系統資訊
docker version > troubleshooting/system_info.txt
docker system info >> troubleshooting/system_info.txt
```

### 效能監控
```bash
# 即時資源監控
docker stats --no-stream

# 磁碟使用分析
du -sh data/ results/ models/ temp/
docker system df

# 記憶體使用分析
free -h
cat /proc/meminfo
```

## 📞 尋求幫助

### 收集診斷資訊
在報告問題前，請收集以下資訊：

```bash
#!/bin/bash
# troubleshooting-info.sh - 收集診斷資訊

echo "=== 系統資訊 ===" > troubleshooting-report.txt
uname -a >> troubleshooting-report.txt
docker version >> troubleshooting-report.txt
docker-compose version >> troubleshooting-report.txt

echo -e "\n=== 容器狀態 ===" >> troubleshooting-report.txt
docker ps -a >> troubleshooting-report.txt

echo -e "\n=== 磁碟空間 ===" >> troubleshooting-report.txt
df -h >> troubleshooting-report.txt
docker system df >> troubleshooting-report.txt

echo -e "\n=== 網路狀態 ===" >> troubleshooting-report.txt
docker network ls >> troubleshooting-report.txt

echo -e "\n=== 錯誤日誌 ===" >> troubleshooting-report.txt
docker-compose -f .devcontainer/docker-compose.yml logs --tail=50 >> troubleshooting-report.txt

echo "診斷資訊已儲存到 troubleshooting-report.txt"
```

### 社群支援
- **GitHub Issues**: 在專案 Repository 建立 Issue
- **文檔更新**: 協助改善此故障排除指南
- **經驗分享**: 分享解決方案供其他使用者參考

## 🎯 預防措施

### 定期維護
```bash
# 每週執行的維護腳本
#!/bin/bash
# weekly-maintenance.sh

echo "🧹 執行每週維護..."

# 1. 清理 Docker 快取
docker system prune -f

# 2. 更新套件
bash .devcontainer/scripts/utils/package-manager.sh update

# 3. 備份環境
bash .devcontainer/scripts/utils/package-manager.sh backup

# 4. 檢查系統狀態
bash .devcontainer/scripts/manage/devcontainer-status.sh

echo "✅ 維護完成"
```

### 最佳實踐
1. **定期備份**: 重要變更前先備份環境
2. **版本鎖定**: 生產環境使用明確的套件版本
3. **資源監控**: 定期檢查磁碟和記憶體使用
4. **日誌保留**: 保留足夠的日誌供故障分析
5. **測試環境**: 在測試環境中驗證變更

---

**遇到問題時保持冷靜，按步驟診斷通常能快速解決**  
*此指南將持續更新以涵蓋更多場景*