# DiHiggs Explorer

A machine-learning–driven toolkit to explore Di-Higgs production in a generic Two-Higgs-Doublet Model (2HDM).  
Efficiently sample a 7-dimensional parameter space, compute Higgs pair branching ratios, filter unphysical points, and build surrogate ML models with interpretability (SHAP).

---

## Table of Contents

- [Features](#features)  
- [Requirements & Installation](#requirements--installation)  
- [Usage](#usage)  
- [Project Structure](#project-structure)  

---

## Features

- **Parameter Sampling**  
  - Latin Hypercube Sampling (LHS), grid or random sampling  
- **Physics Simulation**  
  - Compute Di-Higgs branching ratios (e.g. \(hh \to b\bar b \gamma\gamma\))  
  - Apply theoretical constraints (positivity, unitarity, perturbativity)  
- **Data Cleaning & Analysis**  
  - Identify and filter `NaN` / unphysical points  
  - Exploratory plots: scatter, heatmaps, 3D projections  
- **Machine Learning Models**  
  - Regression (Random Forest, Gaussian Process, Neural Net) for BR prediction  
  - Classification of viable vs. invalid parameter points  
- **Model Interpretability**  
  - SHAP explanations to quantify parameter importance  
- **Active Learning & Optimization**  
  - Bayesian optimization or NSGA-II to propose new sampling points  

---

## Requirements & Installation

### Prerequisites

- **Operating System**: Linux (Debian/Ubuntu), macOS, or Windows (WSL recommended)  
- **Python**: ≥ 3.8  

- **GSL (GNU Scientific Library)**  
  - **Debian/Ubuntu**:  
    ```bash
    dpkg -l | grep libgsl-dev
    ```
    If not installed:
    ```bash
    sudo apt-get update
    sudo apt-get install libgsl-dev
    ```
    Verify header:
    ```bash
    ls /usr/include/gsl/gsl_matrix.h
    ```
  - **macOS (Homebrew)**:
    ```bash
    brew install gsl
    ```
  - **Windows**: use WSL and follow the Debian/Ubuntu instructions above.

### Python Dependencies

All required packages are listed in `requirements.txt`:

- `numpy`, `pandas`, `scipy`  
- `matplotlib`, `seaborn`  
- `scikit-learn`, `xgboost`  
- `shap`  
- `pyDOE`  
- `jupyterlab`  

### Installation Steps

1. **Create and activate a virtual environment**
    
    ```bash
    python3 -m venv venv
    source venv/bin/activate   # Linux/macOS
    venv\Scripts\activate      # Windows
    ```

2. **Install Python dependencies**

    ```bash
    pip install --upgrade pip
    pip install -r requirements.txt
    ```

3. **Recompile C++ libraries**

    ```bash
    # Compile 2HDMC library
    cd 2hdmc/
    make clean
    make
    
    # Compile the main project
    cd ..
    make clean
    make
    ```


## Usage

1. **Parameter Sampling**

   * Open `notebooks/0b1_samplings.ipynb`
   * Configure your sampling strategy (LHS, grid, random)
2. **Simulation & Data Cleaning**

   * Run the sampling notebook to generate raw data
   * Use `scripts/clean_data.py` to filter and label points
3. **Exploratory Analysis**

   * Open `notebooks/1_exploratory_analysis.ipynb` for plots and statistics
4. **Model Training**

   * Open `notebooks/2_ml_models.ipynb`
   * Train regression/classification models, evaluate metrics
5. **Interpretability**

   * Run `notebooks/3_shap_analysis.ipynb` to generate SHAP summary plots
6. **Active Learning Loop**

   * Use `scripts/active_learning.py` to propose new points based on model uncertainty

---

## Project Structure

```
├── notebooks/
│   ├── 0b1_samplings.ipynb
│   ├── 1_exploratory_analysis.ipynb
│   ├── 2_ml_models.ipynb
│   └── 3_shap_analysis.ipynb
├── scripts/
│   ├── clean_data.py
│   └── active_learning.py
├── requirements.txt
├── setup.py
└── README.md
```

