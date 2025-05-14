# ceRNA–scRNA Integration Pipeline

This repository implements an integrated framework for combining competing endogeneous RNA (ceRNA) network data with single-cell RNA sequencing (scRNA-seq) data.  Our dual-view autoencoder leverages graph convolutional and attention mechanism to learn robust latent representaitons that capture both ceRNA regulatory interactions and cell-specific expression profiles.



## 1. Data Preprocessing and ceRNA Network Construction

### ceRNA Network Assembly
The ceRNA network captures regulatory interactions among mRNAs, lncRNAs, and (optionally) circRNAs that compete for shared microRNAs (miRNAs). To build this network:

- **Obtain Interaction Data**  
  - **starBase** – interactions among mRNAs, lncRNAs, and miRNAs (CLIP-Seq)  
  - **LncBase** – lncRNA–miRNA interactions  
  - **miRTarBase** – experimentally validated miRNA–mRNA interactions  
  - **miRcode** – predicted miRNA binding sites on lncRNAs  
  - *(Optional)* circBank / circInteractome – circRNA–miRNA interactions

- **Filtering Interactions**  
  - Apply a confidence score threshold to remove low-confidence edges.

- **Expression-Based Filtering**  
  - Remove any RNA nodes with zero expression across all cells in the scRNA-seq dataset.

### Gene Expression Data Processing
- Start from a normalized scRNA-seq expression matrix.
- Filter to retain only those RNAs (mRNAs, lncRNAs, circRNAs) present in the ceRNA network.
- Result: every network node has an associated single-cell expression profile.

---

## 2. Integrated Encoder Architecture

A dual-view encoder captures both the regulatory network structure and single-cell expression variability.

### Cell-View (scRNA-seq) Branch
1. Build a cell–cell similarity graph (e.g., KNN) from expression profiles.  
2. Apply one or more graph convolution layers to aggregate information among similar cells—this denoises and imputes missing values.

### ceRNA-View Branch
1. Treat the ceRNA network as a graph (nodes = RNAs, edges = regulatory interactions).  
2. Apply graph convolution layers to aggregate regulatory signals across RNAs.

### Graph Attention Mechanism
- Fuse the two views via a graph attention layer.  
- Learn edge-importance weights; use a sigmoid activation (rather than node-wise softmax) to allow global thresholding.  
- Dynamically prune weak/spurious edges, refining the network topology based on reconstruction performance.

---

## 3. Decoder and Loss Formulation

### Dual Decoding Strategy
- **ceRNA Network Reconstruction**  
  Inner-product decoder on the latent RNA embeddings. Use contrastive loss comparing positive (observed) and negative (random) edges.

- **Gene Expression Reconstruction**  
  Multi-layer fully connected decoder to reconstruct the scRNA-seq expression matrix.

### Composite Loss Function
- Combine ceRNA network reconstruction loss and gene expression reconstruction loss.  
- Balance with hyperparameters so neither dominates.  
- Include regularization terms (e.g., L2 penalties) to prevent overfitting.

---

## 4. Graph Pruning and Training

### Adaptive Graph Pruning
- Periodically prune low-quality edges using learned attention coefficients.  
- Ensures the final network reflects true biological interactions.

### Mini-Batch Processing
- For large scRNA-seq datasets, split the cell graph (or cells) into mini-batches.  
- Reduces memory usage while learning a unified embedding space.

### Implementation and Optimization
- **Framework:** PyTorch + torch-geometric  
- **Hardware:** NVIDIA A100 GPUs (or equivalent)  
- **Optimizer:** Adam (tune learning rate, weight decay)  
- **Monitoring:** AUROC on both network and expression reconstruction; adjust hyperparameters via validation.

---

## 5. Databases and Resources Needed

- **ceRNA Interaction Databases**  
  starBase, LncBase, miRTarBase, miRcode  
  *(Optional)* circBank / circInteractome

- **scRNA-seq Data Repositories**  
  GEO (Gene Expression Omnibus)  
  ArrayExpress  
  Single Cell Portal (Broad Institute)

- **Annotation Resources**  
  Ensembl or NCBI RefSeq for gene identifiers and metadata  
  KEGG, Reactome for pathway and functional validation

---

## Summary

This pipeline provides a complete workflow for integrating ceRNA regulatory networks with single-cell RNA-seq data using a dual-view graph neural network. It covers data acquisition and filtering, model architecture (encoder, attention, decoder), adaptive graph pruning, training strategies, and required resources. By jointly modeling post-transcriptional interactions and cell-specific expression, it uncovers cell-type–specific ceRNA mechanisms and enhances our understanding of RNA crosstalk in complex tissues.  
