# RP3
RP3 (RNA Polymerase II Pausing Prediction) is a method for predicting transcriptional elongation Pol II pausing by integrating multi-omics data.
# Motivation
Po II pausing is key regulatory step for implementing precise gene expression programs and holds critical significance in a multitude of essential biological processes. However, the genome-wide profiling of Pol II pausing is currently accessible only in a limited number of cell lines, as experimental measurement techniques remain laborious and costly. Therefore, developing computational prediction methods that leverage the widely available transcriptomic and epigenomic data is desirable strategy for interrogating the information of Pol II pausing for cell types of interest. 
# Method
we develop RP3 (RNA Polymerase II Pausing Prediction), a machine learning framework, to integrate genome sequences, histone modification, gene expression, chromatin accessibility and protein-protein interaction data for predicting Pol II pausing events on the human genome. By taking advantage of Network regularization strategy, RP3 can export a set of biologically meaningful transcription TF cofactors that contribute to elucidating the mechanism underlying Pol pausing. Furthermore, we employed a forward feature selection framework to systematically investigate the combination of histone modification signals, i.e., histone code, that are associated with Pol II pausing.
# Source code
RP3.r: Prediction of Pol II pausing by integrating multi-omics data.
Histone_code.r: Identification of the histone code associated with Pol II pausing
Version 1.0 Last updated: November 29, 2023
# Requirements
R 4.3.2
R package: glmnet, Matrix, foreach, pROC, and PRROC.
# Reference
Lixin Ren, Wanbiao Ma, and Yong Wang, Predicting RNA polymerase II transcriptional elongation pausing and associated histone code. In submission.

