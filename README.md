## 🧬 xEnrich: Organ- and Cell-Level Context Enrichment for Biomarker Discovery

**xEnrich** is an R toolkit for translating differential omics results into biological context.\
It addresses a central question in biomarker discovery:

> *Given a list of biomarkers (differentially expressed genes,proteins), where in the body (organ), in which cell types, and through what systems are these changes most relevant?*

By integrating reference data from the **Human Protein Atlas (IHC)** and **Tabula Sapiens single-cell transcriptomes**, xEnrich maps genes or proteins across tissues and cellular populations to identify enriched organs, dominant cell types, and cross-system patterns.

Designed for translational and systems biology workflows, xEnrich provides lightweight functions and a unified interface that can be used prior to pathway enrichment with tools such as **Reactome**, **enrichR**, or **clusterProfiler**.

The core idea is to bridge statistical discovery and physiological relevance moving away from differential lists to mechanistic hypotheses that may accelerate biomarker validation across complex diseases.
