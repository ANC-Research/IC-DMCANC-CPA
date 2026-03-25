# IC-DMCANC-CPA

**Distributed Multichannel Active Noise Control with Intermittent Communication and Coprocessor Assistance**

---

## 📄 Paper

This repository provides the MATLAB implementation of the method proposed in:

**Junwei Ji**, Dongyuan Shi, Xiaoyi Shen, Woon-Seng Gan, Jie Chen, Jun Yang  
*Implementation of distributed multichannel active noise control with intermittent communication and coprocessor assisted data combination*  
**Mechanical Systems and Signal Processing, 2026**

🔗 https://doi.org/10.1016/j.ymssp.2026.114024

---

## 🚀 Highlights

- 🔹 Distributed multichannel ANC (DMCANC) framework  
- 🔹 Weight-Constrained FxLMS (WCFxLMS) for stability  
- 🔹 Intermittent communication (IC) for practical networks  
- 🔹 Mixed Weight Difference (MWD) for efficient fusion  
- 🔹 Coprocessor-assisted architecture for real-time implementation  

---

## 🧠 Overview

This work proposes an **IC-DMCANC-CPA system**, where:

- Each node performs **local adaptive control** independently  
- Communication occurs only at **intermittent intervals**  
- A **coprocessor** handles global fusion to reduce node burden  

Compared with conventional approaches, the proposed framework:

- ✔ Reduces communication overhead  
- ✔ Maintains system stability under network constraints  
- ✔ Achieves performance close to centralized MCANC  
- ✔ Enables scalable and real-time implementation  

---

## 📂 Repository Structure

```bash
.
├── simulation path/              % Acoustic paths
├── DMANC_CompensateSP.m          % MGDFxLMS (doi:10.1109/TASLPRO.2025.3552932)
├── FedDMCANC_case1.m             % Simulation case 1
├── FedDMCANC_case2.m
├── FedDMCANC_case3.m
├── FedDMCANC_case4.m
├── FedDMCANC_case5.m
├── FedMCANC.m                   % Class and function
├── McANC_FxLMS_SIMO.m           % MEFxLMS
├── compressor_16kHz.mat         % Real recorded noise
└── README.md

```


## ⚙️ Getting Started

### Requirements

- MATLAB (R2021 or newer recommended)
- Signal Processing Toolbox

### ▶️ Run Example

```matlab
run('FedDMCANC_case1.m')
```

## 📊 Results Summary

- Achieves **near-centralized noise reduction performance**
- Significantly reduces **communication frequency**
- Robust under **intermittent and heterogeneous communication**
- Validated through both **simulation and real-time experiments**

---

## 📖 Citation

If you find this work useful, please cite:

```bibtex
@article{ji2026implementation,
title = {Implementation of distributed multichannel active noise control with intermittent communication and coprocessor assisted data combination},
journal = {Mechanical Systems and Signal Processing},
volume = {248},
pages = {114024},
year = {2026},
doi = {https://doi.org/10.1016/j.ymssp.2026.114024},
author = {Junwei Ji and Dongyuan Shi and Xiaoyi Shen and Woon-Seng Gan and Jie Chen and Jun Yang}
}
```

## 📬 Contact

**Junwei Ji**  
Email: JUNWEI002@e.ntu.edu.sg
---

## ⭐ Notes

- This repository focuses on **algorithm validation and simulation**
- Real-time implementation details are described in the paper
- Intended for **research and academic use**
