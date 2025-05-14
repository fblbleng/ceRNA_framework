# ceRNA_framework
# ceRNA–scRNA Integration Framework

An end-to-end computational pipeline and reference implementation for integrating competing endogenous RNA (ceRNA) networks with single-cell RNA-seq (scRNA-seq) data using graph neural networks. This repository provides data-processing scripts, model code, and evaluation tools to uncover cell-type–specific post-transcriptional regulatory interactions.

---

## Features

- **ceRNA Network Assembly**  
  - Harvests high-confidence mRNA–miRNA, lncRNA–miRNA, and optionally circRNA–miRNA interactions from **starBase**, **LncBase**, **miRTarBase**, **miRcode** (and **circBank**/circInteractome).  
  - Filters and harmonizes interaction data into a unified edge list.

- **Data Preprocessing**  
  - Cleans and normalizes scRNA-seq matrices.  
  - Aligns expression profiles to the ceRNA network by removing RNAs with zero expression.

- **Dual-View Graph Encoder**  
  - **Cell-View:** Builds a cell–cell KNN graph to denoise and impute single-cell expression via graph convolutions.  
  - **ceRNA-View:** Applies graph convolutions over the ceRNA network to aggregate regulatory signals.  
  - **Attention Layer:** Learns and prunes low-confidence edges via adaptive, sigmoid-based graph attention.

- **Decoders & Losses**  
  - **Network Reconstruction:** Inner-product decoder with contrastive loss on positive/negative ceRNA edges.  
  - **Expression Reconstruction:** Multi-layer feed-forward decoder optimized with mean squared error on highly variable RNAs.  
  - **Composite Loss:** Balances network and expression reconstruction by tunable hyperparameters.

- **Adaptive Graph Pruning**  
  - Periodically removes edges below a learned global percentile threshold.  
  - Guarantees minimum connectivity per node to avoid disconnection.

- **Scalability & Mini-Batching**  
  - Supports large scRNA-seq datasets by splitting graphs/cells into mini-batches.  
  - Maintains a unified embedding space with minimal batch effects.

- **Evaluation & Benchmarking**  
  - AUROC testing on known functional gene sets (e.g., KEGG pathways) via random-walk propagation.  
  - Topology-free benchmarking against degree-preserving random networks to compute z-scores.

---

## Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/yourusername/ceRNA-scRNA-integration.git
   cd ceRNA-scRNA-integration
