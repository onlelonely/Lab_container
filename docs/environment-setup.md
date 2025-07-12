# 環境設置指南

**版本**: v2.0  
**更新日期**: 2025年7月12日  

## 📋 概覽

本指南詳細說明如何設置和配置 Lab Container 開發環境，包括初始安裝、環境配置、以及常見設置問題的解決方案。

## 🔧 系統需求

### 最低系統需求
- **作業系統**: Ubuntu 20.04+ / CentOS 8+ / macOS 11+ / Windows 10+ (WSL2)
- **記憶體**: 8GB RAM (建議 16GB+)
- **儲存空間**: 50GB 可用空間 (建議 100GB+)
- **Docker**: Docker Engine 20.10+ 或 Docker Desktop
- **網路**: 穩定的網際網路連接

### 建議系統配置
- **處理器**: 8核心 CPU 或以上
- **記憶體**: 32GB RAM 或以上
- **儲存**: SSD 硬碟，200GB+ 可用空間
- **GPU**: NVIDIA GPU (用於 CUDA 加速工具)

## 🐳 Docker 環境設置

### Ubuntu/Debian 安裝 Docker

```bash
# 更新套件列表
sudo apt update

# 安裝必要套件
sudo apt install -y apt-transport-https ca-certificates curl software-properties-common

# 添加 Docker 官方 GPG 金鑰
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

# 添加 Docker 儲存庫
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# 安裝 Docker
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# 啟動 Docker 服務
sudo systemctl start docker
sudo systemctl enable docker

# 將使用者加入 docker 群組
sudo usermod -aG docker $USER
```

## 🚀 Lab Container 安裝

### 克隆儲存庫

```bash
# 克隆專案
git clone https://github.com/onlelonely/Lab_container.git

# 進入專案目錄
cd Lab_container

# 檢查專案結構
ls -la
```

### 環境初始化

```bash
# 賦予執行權限
chmod +x scripts/setup.sh

# 執行初始化腳本
./scripts/setup.sh

# 驗證安裝
./scripts/verify-installation.sh
```

## ⚙️ 環境配置

### Profile 系統配置

```bash
# 建立預設 profile
cp profiles/default.profile.example profiles/my-lab.profile

# 編輯 profile 配置
nano profiles/my-lab.profile
```

## 🔍 驗證安裝

### 系統檢查

```bash
# 檢查 Docker
docker --version

# 檢查 Docker 服務
sudo systemctl status docker

# 測試 Docker 權限
docker run hello-world
```

---

**相關文檔**:
- [Profile 系統說明](./profiles.md)
- [故障排除指南](./troubleshooting.md)