# Efficient Sketch-Based Analysis of Large Single-Cell RNA-seq Datasets

This repository contains code and analysis for comparing various **sketching strategies** used to subsample large-scale single-cell RNA-seq data while preserving biological diversity. We evaluate classical, geometric, and newly proposed methods on the **68k PBMC dataset** from Zheng et al. (2017), focusing on their ability to retain rare cell types, preserve cluster diversity, and execute efficiently.

---

## Project Overview

With the exponential growth of single-cell datasets, sketching methods have become crucial for enabling fast and scalable downstream analysis such as clustering, visualization, and rare cell discovery.

In this project, we:

- Implement and compare **six sketching techniques**:
  - Uniform random sampling
  - Hash-based deterministic sampling
  - Geometric sketching (GeoSketch)
  - $k$-Means medoid sampling
  - Kernel Herding
  - **PCA-guided variance-weighted sampling** (our proposed method)
- Evaluate performance on:
  - **Rare cell capture** (Cluster 12 in 68k PBMCs)
  - **Cluster diversity** retention
  - **Runtime scalability**
- Visualize sketches using **UMAP** and compute quantitative metrics

---

## Dataset

We use the **68k Peripheral Blood Mononuclear Cells (PBMC)** dataset by Zheng et al. (2017) [10x Genomics]. It contains ~68,000 cells profiled using droplet-based RNA-seq.

---

## Repository Structure

.
├── UMaps.ipynb # Main notebook for sketching, plotting and evaluation
├── figures/ # Folder containing generated figures (UMAPs, runtimes, diversity)
├── cluster_diversity_err.png # Cluster coverage barplot
├── time_comparison.png # Runtime comparison for 1000-cell sketch
├── time_scalability.png # Log-scaled runtime vs. sketch size
├── README.md # Project documentation
└── Report.pdf # Full 7-page report with figures, equations, and reference


---

## ⚙️ Sketching Methods & Complexity

| Method         | Description                                                       | Time Complexity                |
|----------------|-------------------------------------------------------------------|--------------------------------|
| Uniform        | Random sampling of cells                                          | $\mathcal{O}(n)$               |
| Hash-Based     | Hash cells by index, sort, select                                 | $\mathcal{O}(N)$               |
| GeoSketch      | Farthest point sampling in feature space                          | $\mathcal{O}(nN)$              |
| $k$-Means      | Cluster and take medoids                                          | $\mathcal{O}(tNk)$             |
| Kernel Herding | Greedy selection to match mean embedding                          | $\mathcal{O}(nN)$              |
| PCA-Guided     | Sample based on PCA vector norms                                  | $\mathcal{O}(Nd + n)$          |

Where:
- $N$ = number of cells (68,000)
- $n$ = sketch size
- $d$ = number of PCA dimensions
- $t$ = number of iterations for clustering

---

## 🧠 PCA-Guided Sampling (New Method)

We propose a **simple yet effective** method that leverages PCA scores to preferentially sample cells with higher variance, improving the chances of capturing rare transcriptional states while maintaining speed. Sampling is done proportionally to squared Euclidean norm of PCA vectors.

### Implementation:
```python
variances = np.sum(X_pca**2, axis=1)
probabilities = variances / np.sum(variances)
pca_indices = np.random.choice(len(X_pca), size=n, replace=False, p=probabilities)





