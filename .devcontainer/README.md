# Dev Container Configuration

This directory contains the complete development container configuration for data science projects using Python, R, and bioinformatics tools.

## Directory Structure

```
.devcontainer/
├── devcontainer.json          # VS Code dev container configuration
├── Dockerfile                 # Container build instructions
├── docker-compose.yml         # Multi-container orchestration
├── .env.template             # Environment variables template
├── README.md                 # This file
├── configs/                  # Configuration files
│   ├── python/
│   │   ├── requirements-core.txt
│   │   ├── requirements-bioinformatics.txt
│   │   └── requirements-extended.txt
│   ├── r/
│   │   ├── core-packages.R
│   │   ├── bioconductor-packages.R
│   │   └── extended-packages.R
│   └── conda/
│       ├── environment-core.yml
│       └── environment-extended.yml
├── scripts/                  # Installation and management scripts
│   ├── install/
│   │   ├── install-core.sh
│   │   └── install-extended.sh
│   ├── manage/               # Package management scripts (future)
│   └── utils/                # Utility functions
│       ├── logging.sh
│       ├── retry-logic.sh
│       └── validation.sh
└── tests/                    # Test suites
    ├── test-installation.sh
    └── test-environment.sh
```

## Features

### Security
- Non-root user (devuser) with sudo access
- Proper file permissions and ownership
- Secure container configuration

### Package Management
- Modular configuration files for different package types
- Retry logic with exponential backoff
- Comprehensive error handling and logging
- Validation of installations

### Development Tools
- Python 3.11 with data science libraries
- R 4.3 with tidyverse and Bioconductor
- Jupyter Lab with Python and R kernels
- VS Code extensions for Python, R, and Jupyter

### Testing
- Automated installation testing
- Environment validation
- Performance benchmarking
- Comprehensive test reports

## Usage

### First Time Setup

1. Copy environment template:
   ```bash
   cp .devcontainer/.env.template .devcontainer/.env
   ```

2. Customize environment variables in `.env` file

3. Open in VS Code with Dev Containers extension

### Running Tests

```bash
# Test installation
bash .devcontainer/tests/test-installation.sh

# Test environment
bash .devcontainer/tests/test-environment.sh
```

### Package Management

Core packages are installed automatically. Extended packages can be controlled via environment variables:

```bash
export INSTALL_BIOCONDUCTOR=true
export INSTALL_EXTENDED_PACKAGES=false
```

### Logs

All logs are stored in `/workspace/logs/` with timestamped filenames:
- Installation logs: `devcontainer-YYYYMMDD.log`
- Test results: `test-results/`

## Configuration Files

### Core Dependencies
- `configs/conda/environment-core.yml`: Essential conda packages
- `configs/python/requirements-core.txt`: Core Python packages
- `configs/r/core-packages.R`: Essential R packages

### Extended Dependencies
- `configs/conda/environment-extended.yml`: Additional conda packages
- `configs/python/requirements-extended.txt`: Extended Python packages
- `configs/python/requirements-bioinformatics.txt`: Bioinformatics tools
- `configs/r/bioconductor-packages.R`: Bioconductor packages
- `configs/r/extended-packages.R`: Additional R packages

## Environment Variables

Key environment variables (see `.env.template` for full list):

- `PYTHON_VERSION`: Python version (default: 3.11)
- `R_VERSION`: R version (default: 4.3)
- `INSTALL_BIOCONDUCTOR`: Install Bioconductor packages (default: true)
- `INSTALL_EXTENDED_PACKAGES`: Install extended packages (default: false)
- `PARALLEL_JOBS`: Number of parallel installation jobs (default: 4)

## Troubleshooting

### Common Issues

1. **Package installation failures**: Check logs in `/workspace/logs/`
2. **Network connectivity**: Verify internet connection and proxy settings
3. **Permission issues**: Ensure proper user permissions are set
4. **Resource constraints**: Check available memory and disk space

### Test Commands

```bash
# Quick validation
bash .devcontainer/scripts/utils/validation.sh

# Full test suite
bash .devcontainer/tests/test-installation.sh
bash .devcontainer/tests/test-environment.sh
```

## Development

### Adding New Packages

1. Add packages to appropriate configuration files
2. Update installation scripts if needed
3. Add validation tests
4. Update documentation

### Modifying Scripts

All scripts use the common utility functions:
- `logging.sh`: Standardized logging
- `retry-logic.sh`: Retry mechanisms
- `validation.sh`: Testing and validation

## Performance

Expected installation times:
- Core packages: 5-10 minutes
- Extended packages: 15-30 minutes
- Bioconductor: 20-40 minutes

Resource requirements:
- Memory: 4GB+ recommended
- Disk space: 10GB+ recommended
- CPU: 2+ cores recommended