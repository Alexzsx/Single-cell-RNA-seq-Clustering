# Review for single-cell RNA-seq clustering
Recently, the advanced development of single-cell RNA-seq techniques has become a powerful tool to measure large-scale transcriptomic profiling in a single-cell resolution. Unsupervised clustering has become the central component of single-cell RNA-seq interpretative approach to identifying rare or new cell types, gene expression patterns, and functional implications of stochastic transcription. In this study, we review several popular applications of clustering methods including {\it k}-means clustering, hierarchical clustering, community-detection-based clustering, and density-based clustering methods for its advantages and limitations on single-cell RNA-seq data. In addition, we also review the upstream steps of single-cell RNA-seq clustering such as quality control, normalization, and dimension reduction that have a substantial impact to speed up further calculations. We then conducted several experiments to evaluate and compare the performance and efficiency of several popular single-cell RNA-seq clustering approaches on two transcriptomic profiling datasets.

# PREREQUISITE
It was conducted by R v3.5.3. 

Following R packages should be installed:
<ul>
<li><p>cidr</p></li>
<li><p>SC3</p></li>
<li><p>Seurat</p></li>
<li><p>RaceID3</p></li>
<li><p>SIMLR</p></li>
<li><p>monocle</p></li>  
</ul>

# CONTAINS:
<ul>
<li><p>Clusstering for Human Brain.R : R script to evaluate clustering methods on Human Brain scRNA-Seq Dataset</p></li>
<li><p>Clusstering for Clusstering for.R : R script to evaluate clustering methods on Human Pancreatic Islet scRNA-Seq Dataset</p></li> 
</ul>

## Human Brain scRNA-Seq Dataset
In this dataset there are 420 cells in 8 cell types after we exclude hybrid cells.

Reference for the human brain dataset:

Darmanis, S. et al. A survey of human brain transcriptome diversity at the single cell level. Proceedings of the National Academy of Sciences 112, 7285–7290 (2015).

## Human Pancreatic Islet scRNA-Seq Dataset
In this dataset there are 60 cells in 6 cell types after we exclude undefined cells and bulk RNA-Seq samples.

Reference for the human pancreatic islet dataset:

Li, J. et al. Single-cell transcriptomes reveal characteristic features of human pancreatic islet cell types. EMBO Reports 17, 178–187 (2016).
