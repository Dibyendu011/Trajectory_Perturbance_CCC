# Single-Cell Multi-Omics Analysis Toolkit

This repository contains three independent computational pipelines for advanced single-cell RNA sequencing analysis, covering cell-cell communication, gene regulatory network perturbation, and trajectory inference.

## Analysis Modules

### 1. Cell-Cell Communication Analysis
**File:** `cell_cell_communication.ipynb`  
**Method:** LIANA (LIgand-receptor ANalysis frAmework)  
**Input Data:** `kang.h5ad` - Single-cell RNA-seq dataset  
**Description:** Identifies and quantifies ligand-receptor interactions between different cell types using multiple consensus methods.

### 2. Gene Regulatory Network Perturbation Analysis
**File:** `perturbation.ipynb`  
**Method:** CellOracle  
**Input Data:** 
- **scRNA-seq:** scVelo pancreas development dataset
- **scATAC-seq:** Mouse single-cell ATLAS base GRN (`celloracle.data.load_mouse_scATAC_atlas_base_GRN()`)

**Key Objects:**
- `links_object.celloracle.links` - Contains Gene Regulatory Network connections (TF → target gene relationships)
- `my_oracle.celloracle.oracle` - Main analysis object storing GRN model and perturbation simulations

**Description:** Models transcription factor perturbations and predicts downstream effects on cell identity and state transitions.

### 3. Trajectory Inference Analysis
**File:** `trajectory.R`  
**Method:** Monocle3  
**Input Data:** [GSE122662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122662) - Single-cell RNA sequencing of developing pancreas  
**Description:** Reconstructs developmental trajectories and pseudotemporal ordering of cells during pancreas development.

## Dataset Sources

- **Cell-Cell Communication:** Processed Kang et al. dataset (`kang.h5ad`)
- **GRN Perturbation:** 
  - scRNA-seq: Pancreas development data from scVelo
  - scATAC-seq: Mouse chromatin accessibility atlas
- **Trajectory Inference:** Reconstructs developmental trajectories and pseudotemporal ordering of cells from [GSE122662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122662) during cellular reprogramming.
