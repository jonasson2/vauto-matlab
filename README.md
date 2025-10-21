# vauto – Exact VARMA Likelihood in MATLAB

**vauto** is a MATLAB library for computing the exact likelihood and gradient of VARMA (Vector Autoregressive Moving Average) models, supporting both complete and incomplete data.  
It accompanies **Algorithm 878** and its companion methodology paper, published in *ACM Transactions on Mathematical Software (TOMS), 35(1), 2008*.

This repository hosts the actively maintained MATLAB version (currently **v1.1.0**).  
Earlier releases (**v1.0.0**, **v1.0.1**) are available under **GitHub Releases**.

---

## ✨ Features

- Exact likelihood and gradient evaluation for VARMA models  
- Support for missing data  
- Numerical stability via block matrix factorisation  
- Unit-test suite (`test/`) for verification  
- Includes required optimization helpers (`immoptibox/`)

---

## 📦 Installation

1. Clone or download the repository:
   ```bash
   git clone https://github.com/jonasson2/vauto-matlab.git